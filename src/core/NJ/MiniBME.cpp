#include "MiniBME.hpp"
#include <NJ/MiniNJ.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <search/UNNISearch.hpp>
#include <limits>

using DistanceVectorMatrix = std::vector<DistanceMatrix>;

static BoolMatrix getBoolMatrix(unsigned int N,
    unsigned int M,
    bool value = false)
{
  std::vector<bool> v(M, value);
  return BoolMatrix(N, v);
}

static DistanceMatrix getNullMatrix(unsigned int N,
    double value = 0.0)
{
  std::vector<double> nullDistances(N, value);
  return DistanceMatrix(N, nullDistances); 
}
static DistanceMatrix getMatrix(unsigned int N,
    unsigned int M,
    double value)
{
  std::vector<double> nullDistances(M, value);
  return DistanceMatrix(N, nullDistances); 
}


static void fillDistancesRec(pll_unode_t *currentNode, 
    const std::vector<unsigned int> &nodeIndexToSpid,
    double currentDistance,
    std::vector<double> &distances,
    std::vector<bool> *belongsToPruned)
{
  if (!belongsToPruned || (*belongsToPruned)[currentNode->node_index]) {
    currentDistance += 1.0;
  }
  if (!currentNode->next) {
    // leaf
    if (!belongsToPruned || (*belongsToPruned)[currentNode->node_index]) {
      distances[nodeIndexToSpid[currentNode->node_index]] 
        = currentDistance;
    }
    return;
  }
  fillDistancesRec(currentNode->next->back, 
      nodeIndexToSpid,
      currentDistance, 
      distances,
      belongsToPruned);
  fillDistancesRec(currentNode->next->next->back, 
      nodeIndexToSpid,
      currentDistance, 
      distances,
      belongsToPruned);
}

static void fillSpeciesDistances(const PLLUnrootedTree &speciesTree,
    const StringToUint &speciesStringToSpeciesId,
    DistanceMatrix &speciesDistanceMatrix,
    std::vector<bool> *belongsToPruned = nullptr)
{
  std::vector<unsigned int> nodeIndexToSpid(speciesTree.getLeavesNumber());
  for (auto leaf: speciesTree.getLeaves()) {
    nodeIndexToSpid[leaf->node_index] = speciesStringToSpeciesId.at(leaf->label);
  }
  for (auto leaf: speciesTree.getLeaves()) {
    if (belongsToPruned && !(*belongsToPruned)[leaf->node_index]) {
      continue;
    }
    std::string label = leaf->label;
    auto id = speciesStringToSpeciesId.at(label);
    fillDistancesRec(leaf->back,
        nodeIndexToSpid,
        0.0,
        speciesDistanceMatrix[id],
        belongsToPruned);
  }
}

static void getPrunedSpeciesMatrix(const PLLUnrootedTree &speciesTree,
  const StringToUint &speciesStringToSpeciesId,
  const std::unordered_set<std::string> &coverage,
  DistanceMatrix &distanceMatrix)
{
  std::vector<bool> hasChildren(speciesTree.getDirectedNodesNumber(), false);
  std::vector<bool> belongsToPruned(speciesTree.getDirectedNodesNumber(), false);
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto index = node->node_index;
    if (node->next) {
      auto leftIndex = node->next->back->node_index;
      auto rightIndex = node->next->next->back->node_index;
      hasChildren[index] = hasChildren[leftIndex] || hasChildren[rightIndex];
      belongsToPruned[index] = hasChildren[leftIndex] && hasChildren[rightIndex];
    } else {
      bool isCovered = coverage.find(std::string(node->label)) != coverage.end();
      hasChildren[index] = isCovered;
      belongsToPruned[index] = isCovered;
    }
  }
  fillSpeciesDistances(speciesTree,
      speciesStringToSpeciesId,
      distanceMatrix,
      &belongsToPruned);
}
    




MiniBME::MiniBME(const PLLUnrootedTree &speciesTree,
    const Families &families,
    bool pruneMode):
  _prune(pruneMode)
{
  bool minMode = true;
  bool reweight = false;
  bool ustar = false;
  for (auto leaf: speciesTree.getLeaves()) {
    std::string label(leaf->label);
    _speciesStringToSpeciesId.insert({label, _speciesIdToSpeciesString.size()});
    _speciesIdToSpeciesString.push_back(label);
  }
  if (_prune) {
    // fill species-spid mappings
    PerCoreGeneTrees::getPerCoreFamilies(families, _perCoreFamilies);
    _geneDistanceMatrices.resize(_perCoreFamilies.size());
    _geneDistanceDenominators.resize(_perCoreFamilies.size());
    unsigned int speciesNumber = speciesTree.getLeavesNumber();
    _perFamilyCoverageStr.resize(_perCoreFamilies.size());
    _perFamilyCoverage.resize(_perCoreFamilies.size());
    // map each species string to a species ID
    for (unsigned int i = 0; i < _perCoreFamilies.size(); ++i) {
      auto &family = _perCoreFamilies[i];
      _perFamilyCoverage[i] = std::vector<bool>(speciesTree.getLeavesNumber(), false);
      _geneDistanceMatrices[i] = getNullMatrix(speciesNumber);
      _geneDistanceDenominators[i] = getNullMatrix(speciesNumber);
      GeneSpeciesMapping mappings;
      mappings.fill(family.mappingFile, family.startingGeneTree);
      for (auto species: mappings.getCoveredSpecies()) {
        _perFamilyCoverageStr[i].insert(species);
        _perFamilyCoverage[i][_speciesStringToSpeciesId.at(species)] = true;
      }
      std::ifstream reader(family.startingGeneTree);
      std::string geneTreeStr;
      while (std::getline(reader, geneTreeStr)) {
        PLLUnrootedTree geneTree(geneTreeStr, false);
        MiniNJ::geneDistancesFromGeneTree(geneTree, 
            mappings,
            _speciesStringToSpeciesId,
            _geneDistanceMatrices[i],
            _geneDistanceDenominators[i],
            minMode,
            reweight,
            ustar);
      }
    }
    _prunedSpeciesMatrices = std::vector<DistanceMatrix>(
        _perCoreFamilies.size(), getNullMatrix(speciesNumber));
  } else {  
    _geneDistanceMatrices.resize(1);
    MiniNJ::computeDistanceMatrix(families,
      minMode, 
      reweight,
      ustar,
      _geneDistanceMatrices[0],
      _speciesIdToSpeciesString,
      _speciesStringToSpeciesId);
  }
}


double MiniBME::computeBME(const PLLUnrootedTree &speciesTree)
{
  if (_prune) {
    return _computeBMEPrune(speciesTree);
  }
  unsigned int N = _speciesIdToSpeciesString.size();
  DistanceMatrix speciesDistanceMatrix = getNullMatrix(N);
  fillSpeciesDistances(speciesTree, 
      _speciesStringToSpeciesId,
      speciesDistanceMatrix);
  double res = 0.0;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      double weight = pow(2.0, speciesDistanceMatrix[i][j]);
      res += _geneDistanceMatrices[0][i][j] / weight;
    }
  }
  _prunedSpeciesMatrices.resize(1);
  _prunedSpeciesMatrices[0] = speciesDistanceMatrix;
  _computeSubBMEs(speciesTree);
  return res;
}


double MiniBME::_computeBMEPrune(const PLLUnrootedTree &speciesTree)
{
  unsigned int N = _speciesIdToSpeciesString.size();
  DistanceMatrix speciesDistanceMatrix = getNullMatrix(N);
  fillSpeciesDistances(speciesTree, 
      _speciesStringToSpeciesId,
      speciesDistanceMatrix);
  double res = 0.0;
  for (unsigned int k = 0; k < _perCoreFamilies.size(); ++k) {
    getPrunedSpeciesMatrix(speciesTree, 
          _speciesStringToSpeciesId,
          _perFamilyCoverageStr[k],
          _prunedSpeciesMatrices[k]);    
  }
  for (unsigned k = 0; k < _geneDistanceMatrices.size(); ++k) {
    for (unsigned int i = 0; i < N; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        if (0.0 == _geneDistanceDenominators[k][i][j]) {
          continue;
        }
        res += _geneDistanceMatrices[k][i][j] / pow(2.0, _prunedSpeciesMatrices[k][i][j]);
      }
    }
  }
  ParallelContext::sumDouble(res);
  _computeSubBMEsPrune(speciesTree);
  return res;
}



bool isNumber(double v) {
  return v == 0.0 || std::isnormal(v);
}


void MiniBME::_computeSubBMEs(const PLLUnrootedTree &speciesTree)
{
  auto subtrees1 = speciesTree.getReverseDepthNodes();
  auto nodesNumber = subtrees1.size();
  _subBMEs = DistanceVectorMatrix(nodesNumber, 
      getMatrix(nodesNumber, 1, std::numeric_limits<double>::infinity()));
  for (auto n1: subtrees1) {
    auto subtrees2 = speciesTree.getPostOrderNodesFrom(n1->back);
    auto i1 = n1->node_index;
    for (auto n2: subtrees2) {
      auto i2 = n2->node_index;
      if (!n1->next && !n2->next) { // both leaves
        // edge case
        _subBMEs[i1][i2][0] = _geneDistanceMatrices[0][i1][i2];
      } else if (n1->next && !n2->next) { // n2 is a leaf
        // we already computed the symetric
        _subBMEs[i1][i2][0] = _subBMEs[i2][i1][0];
      } else  { // n2 is not a leaf
        auto left2 = n2->next->back->node_index;
        auto right2 = n2->next->next->back->node_index;
        _subBMEs[i1][i2][0] = 0.5 * (_subBMEs[i1][left2][0] + _subBMEs[i1][right2][0]);
      }
    }
  }

}


// computes _subBMEs[0][i1][i2]
static void computeSubBMEsPruneRec(pll_unode_t *n1,
    pll_unode_t *n2,
    BoolMatrix &treated,
    BoolMatrix &hasChildren,
    BoolMatrix &belongsToPruned,
    std::vector<DistanceMatrix> &geneDistances,
    DistanceVectorMatrix &subBMEs)
{
  auto i1 = n1->node_index;
  auto i2 = n2->node_index;
  unsigned int K = subBMEs[0][0].size();
  if (treated[i1][i2]) {
    // stop if we've already been there
    return;
  }
  if (!n1->next && !n2->next) {
    // leaf-leaf case
    for (unsigned int k = 0; k < K; ++k) {
      if (belongsToPruned[i1][k] && belongsToPruned[i2][k]) {
        subBMEs[i1][i2][k] = geneDistances[k][i1][i2];
      } else {
        subBMEs[i1][i2][k] = 0.0;
      }
    }
  } else if (!n2->next) {
    // only n2 is a leaf, compute the symetric
    computeSubBMEsPruneRec(n2, n1, treated, hasChildren, 
        belongsToPruned, geneDistances, subBMEs);
    subBMEs[i1][i2] = subBMEs[i2][i1];
  } else {
    // n2 is not a leaf, we can run recursion on n2
    auto left2 = n2->next->back;
    auto right2 = n2->next->next->back;
    computeSubBMEsPruneRec(n1, left2, treated, hasChildren, 
        belongsToPruned, geneDistances, subBMEs);
    computeSubBMEsPruneRec(n1, right2, treated, hasChildren, 
        belongsToPruned, geneDistances, subBMEs);
    for (unsigned int k = 0; k < K; ++k) {
      if (!hasChildren[i1][k] || !hasChildren[i2][k]) {
        // Edge case: either n1 or n2 is not in the path
        // of the induced tree, so its subBME should be
        // ignored
        subBMEs[i1][i2][k] = 0.0;
        continue;
      }
      // Now we can assume that both n1 and n2 are in 
      // the path of the induced tree, but they might not
      // be binary nodes in the induced tree
      if (!belongsToPruned[i2][k]) {
        // Case 1: n2 is not in the induced tree
        // Go down to n2's child that has children
        if (hasChildren[left2->node_index][k]) {
          subBMEs[i1][i2][k] = subBMEs[i1][left2->node_index][k];
          assert(isNumber(subBMEs[i1][i2][k]));
        } else {
          assert(hasChildren[right2->node_index][k]);
          subBMEs[i1][i2][k] = subBMEs[i1][right2->node_index][k];
          assert(isNumber(subBMEs[i1][i2][k]));
        }
      } else if (!belongsToPruned[i1][k]) {
        // Case 2: n2 is in the induced tree, n1 is not
        // Note that n1 cannot be a leaf at this point
        // go down to n1's child that has children
        assert(n1->next != nullptr);
        auto left1 = n1->next->back;
        auto right1 = n1->next->next->back;
        auto lefti1 = n1->next->back->node_index;
        auto righti1 = n1->next->next->back->node_index;
        computeSubBMEsPruneRec(left1, n2, treated, hasChildren, 
            belongsToPruned, geneDistances, subBMEs);
        computeSubBMEsPruneRec(right1, n2, treated, hasChildren, 
            belongsToPruned, geneDistances, subBMEs);
        if (hasChildren[lefti1][k]) {
          subBMEs[i1][i2][k] = subBMEs[lefti1][i2][k];
          assert(isNumber(subBMEs[i1][i2][k]));
        } else {
          assert(hasChildren[righti1][k]);
          subBMEs[i1][i2][k] = subBMEs[righti1][i2][k];
          assert(isNumber(subBMEs[i1][i2][k]));
        }
      } else {
        // case 3: both n1 and n2 are in the induced tree
        subBMEs[i1][i2][k] = 0.5 * (
            subBMEs[i1][left2->node_index][k] + 
            subBMEs[i1][right2->node_index][k]);
        assert(isNumber(subBMEs[i1][i2][k]));
      }
  
    }
  }


  for (unsigned k = 0; k < K; ++k) {
    assert(isNumber(subBMEs[i1][i2][k]));
  }
  treated[i1][i2] = true;
}

void MiniBME::_computeSubBMEsPrune(const PLLUnrootedTree &speciesTree)
{
  auto nodesNumber = speciesTree.getDirectedNodesNumber();
  auto K = _perCoreFamilies.size();
  // Fill hasChildren and belongsToPruned
  _hasChildren = getBoolMatrix(nodesNumber, 
      _perCoreFamilies.size(),
      false);
  _belongsToPruned = _hasChildren;
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto index = node->node_index;
    for (unsigned int k = 0; k < K; ++k) {
      if (!node->next) {
        if (_perFamilyCoverage[k][index]) {
          _hasChildren[index][k] = true;
          _belongsToPruned[index][k] = true;
        }
      } else {
        auto left = node->next->back->node_index;
        auto right = node->next->next->back->node_index;
        _hasChildren[index][k] = 
          _hasChildren[left][k] || _hasChildren[right][k];
        _belongsToPruned[index][k] = 
          _hasChildren[left][k] && _hasChildren[right][k];
      }
    }
  }
  // Fill  the per-family subBME matrices
  BoolMatrix treated = getBoolMatrix(nodesNumber, nodesNumber, false);
  _subBMEs = DistanceVectorMatrix(nodesNumber, 
      getMatrix(nodesNumber, 
        K, 
        std::numeric_limits<double>::infinity()));
  for (auto n1: speciesTree.getPostOrderNodes()) {
    
    // we only need the subBMEs of the nodes of the
    // subtree rooted at n1->back
    for (auto n2: speciesTree.getPostOrderNodesFrom(n1->back)) {
      if (n2 == n1->back) {
        continue;
      }
      computeSubBMEsPruneRec(n1,
        n2,
        treated,
        _hasChildren,
        _belongsToPruned,
        _geneDistanceMatrices,
        _subBMEs);
    }
  }
}

double MiniBME::computeNNIDiff(const PLLUnrootedTree &speciesTree,
      const UNNIMove &nni)
{
  auto A = nni.getA()->node_index;
  auto B = nni.getB()->node_index;
  auto C = nni.getC()->node_index;
  auto D = nni.getD()->node_index;
  auto e1 = nni.edge->node_index;;
  auto e2 = nni.edge->back->node_index;;
  if (!_prune) {
    auto diffPlus = _subBMEs[A][C][0] + _subBMEs[B][D][0];
    auto diffMinus = _subBMEs[A][B][0] + _subBMEs[C][D][0];
    return (diffPlus - diffMinus) * 0.125;
  }

  double diffPlus = 0.0;
  double diffMinus = 0.0;
  for (unsigned int k = 0; k < _perCoreFamilies.size(); ++k) {
    if (!_belongsToPruned[e1][k] || !_belongsToPruned[e2][k]) {  
    // skip families for which the NNI move does not change
      // the induced tree topology, because the formula
      // does not hold anymore in this case
      continue;
    }
    diffPlus += _subBMEs[A][C][k] + _subBMEs[B][D][k];
    diffMinus += _subBMEs[A][B][k] + _subBMEs[C][D][k];
  }
  double res = (diffPlus - diffMinus) * 0.125;
  ParallelContext::sumDouble(res);
  return res;
}
