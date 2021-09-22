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

// return the maximum non-diagonal value, assuming
// that m is symetric
double getMaxSym(const DistanceMatrix &m) {
  double res = -std::numeric_limits<double>::infinity();
  for (unsigned int i = 0; i < m.size(); ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      res = std::max(res, m[i][j]);
      Logger::info << m[i][j] << std::endl;
    }
  }
  return res;
}

double MiniBME::_computeBMEPrune(const PLLUnrootedTree &speciesTree)
{
  unsigned int N = _speciesIdToSpeciesString.size();
  DistanceMatrix speciesDistanceMatrix = getNullMatrix(N);
  fillSpeciesDistances(speciesTree, 
      _speciesStringToSpeciesId,
      speciesDistanceMatrix);
  //double maxSpeciesDistance = getMaxSym(speciesDistanceMatrix);
  double res = 0.0;
  // O(kn^2)
  for (unsigned int k = 0; k < _perCoreFamilies.size(); ++k) {
    // O(n^2)
    getPrunedSpeciesMatrix(speciesTree, 
          _speciesStringToSpeciesId,
          _perFamilyCoverageStr[k],
          _prunedSpeciesMatrices[k]);    
  }
  // O(kn^2)
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

double MiniBME::computeNNIDiff(const UNNIMove &nni)
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

static pll_unode_t *getOtherNext(pll_unode_t *n1, 
    pll_unode_t *n2)
{
  if (n1->next == n2) {
    return n2->next;
  } else {
    assert(n1->next->next == n2);
    return n1->next;
  }
}

// test regrafting  Wp between Vs and Vsplus1 after 
// having regrafted Wp between Vsminus and Vs
// and recursively test further SPR mvoes
static void getBestSPRRec(unsigned int s,
    pll_unode_t *W0, 
    pll_unode_t *Wp, 
    pll_unode_t *Wsminus1, 
    pll_unode_t *Vsminus1, 
    double delta_Vsminus2_Wp, // previous deltaAB
    pll_unode_t *Vs, 
    double Lsminus1, // L_s-1
    pll_unode_t *&bestRegraftNode,
    double &bestLs,
    const std::vector<DistanceMatrix> &subBMEs)
{
  pll_unode_t *Ws = getOtherNext(Vs, Vsminus1->back)->back; 
  // compute L_s
  // we compute LS for an NNI to ((A,B),(C,D)) by swapping B and C
  // where:
  // - A is Vsminus1
  // - B is Wp
  // - C is Ws
  // - D is Vs->back
  // A was modified by the previous (simulated) NNI moves
  // and its average distances need an update
  double diff = 0.0;
  double deltaCD = subBMEs[Ws->node_index][Vs->back->node_index][0];
  double deltaBD = subBMEs[Wp->node_index][Vs->back->node_index][0];
  double deltaAB = 0.0;
  double deltaAC = 0.0;
  if (s == 1) {
    deltaAB = subBMEs[W0->node_index][Wp->node_index][0];
    deltaAC = subBMEs[W0->node_index][Ws->node_index][0];
  } else {
    deltaAB = 0.5 * (delta_Vsminus2_Wp + 
      subBMEs[Wsminus1->node_index][Wp->node_index][0]);
    deltaAC = subBMEs[Vsminus1->node_index][Ws->node_index][0];
    deltaAC -= pow(0.5, s) * subBMEs[Wp->node_index][Ws->node_index][0];
    deltaAC += pow(0.5, s) * subBMEs[W0->node_index][Ws->node_index][0];
  }
  diff = 0.125 * (deltaAB + deltaCD - deltaAC - deltaBD);
  double Ls = Lsminus1 + diff;
  if (Ls > bestLs) {
    bestLs = Ls;
    bestRegraftNode = Vs;
  }
  // recursive call
  if (Vs->back->next) {
    getBestSPRRec(s+1, W0, Wp, Ws, Vs, deltaAB, Vs->back->next, 
        Ls, bestRegraftNode, bestLs, subBMEs);
    getBestSPRRec(s+1, W0, Wp, Ws, Vs, deltaAB, Vs->back->next->next, 
        Ls, bestRegraftNode, bestLs, subBMEs);
  }
}


static void getBestSPRRecMissing(unsigned int s,
    std::vector<unsigned int> sprime, // copy!!
    std::vector<pll_unode_t *> W0s, 
    pll_unode_t *Wp, 
    pll_unode_t *Wsminus1, 
    pll_unode_t *Vsminus1, 
    std::vector<double> delta_Vsminus2_Wp, // previous deltaAB
    pll_unode_t *Vs, 
    double Lsminus1, // L_s-1
    pll_unode_t *&bestRegraftNode,
    double &bestLs,
    unsigned int &bestS,
    const std::vector<DistanceMatrix> &subBMEs,
    const BoolMatrix &belongsToPruned,
    const BoolMatrix &hasChildren,
    std::vector<bool> Vsminus2HasChildren, // does Vsminus2 have children after the previous moves
    std::vector<bool> Vsminus1HasChildren) // does Vsminus1 have children after the previous moves
{
  std::vector<double> pows(s+1);
  for (unsigned int i = 0; i <= s; ++i) {
    pows[i] = pow(0.5, i);
  }
  unsigned int maxRadius = 9999;
  pll_unode_t *Ws = getOtherNext(Vs, Vsminus1->back)->back; 
  unsigned int K = subBMEs[0][0].size();
  // compute L_s
  // we compute LS for an NNI to ((A,B),(C,D)) by swapping B and C
  // where:
  // - A is Vsminus1
  // - B is Wp
  // - C is Ws
  // - D is Vs->back
  // A was modified by the previous (simulated) NNI moves
  // and its average distances need an update
  double diff = 0.0;
  std::vector<double> delta_Vsminus1_Wp(K, std::numeric_limits<double>::infinity());
  std::vector<bool> VsHasChildren = Vsminus1HasChildren;
  for (unsigned int k = 0; k < K; ++k) {
    // if Wp has no children in the induced tree, the SPR move has no effect
    // on the score of this family
    if (!hasChildren[Wp->node_index][k]) {
      continue;
    }
    // if there is no regrafting branch in the induced tree anymore, there is no 
    // point in continuing computing stuff for this family
    if (!hasChildren[Vsminus1->back->node_index][k]) {
      continue;
    }
    // The current NNI move affects the current family score if and only if each of
    // its 4 subtrees (after having recursively applied the previous NNI moves!!)
    // have children in the induced tree
    // At this point, Wp has children. Vsminus1HasChildren tells if the updated Vs has children
    // belongsToPruned == true iif both the node children (Ws and Vs->back) have children
    bool applyDiff = Vsminus1HasChildren[k] && belongsToPruned[Vsminus1->back->node_index][k];
 
    // only matters from s == 3
    if (nullptr == W0s[k] && Vsminus1HasChildren[k]) {
      W0s[k] = Wsminus1;
    }

    // deltaCD and deltaBD are trivial because not affected by previous NNI moves
    double deltaCD = subBMEs[Ws->node_index][Vs->back->node_index][k];
    double deltaBD = subBMEs[Wp->node_index][Vs->back->node_index][k];
    
    // deltaAB and deltaAC are affected by the previous moves
    double deltaAB = 0.0;
    double deltaAC = 0.0;
    if (sprime[k] <= 1 && applyDiff) {
      deltaAB = subBMEs[W0s[k]->node_index][Wp->node_index][k];
      deltaAC = subBMEs[W0s[k]->node_index][Ws->node_index][k];
    } else if (W0s[k]) {
      deltaAC = subBMEs[Vsminus1->node_index][Ws->node_index][k];
      deltaAC -= pows[sprime[k]] * subBMEs[Wp->node_index][Ws->node_index][k];
      deltaAC += pows[sprime[k]] * subBMEs[W0s[k]->node_index][Ws->node_index][k];
      if (hasChildren[Wsminus1->node_index][k] && Vsminus2HasChildren[k]) {
        deltaAB = 0.5 * (delta_Vsminus2_Wp[k] + 
          subBMEs[Wsminus1->node_index][Wp->node_index][k]);
      } else if (hasChildren[Wsminus1->node_index][k] && !Vsminus2HasChildren[k]) {
        deltaAB = subBMEs[Wsminus1->node_index][Wp->node_index][k];
      } else if (!hasChildren[Wsminus1->node_index][k] && Vsminus2HasChildren[k]) {
        deltaAB = delta_Vsminus2_Wp[k];
      } else {
        assert(!applyDiff);
        deltaAB = std::numeric_limits<double>::infinity(); // won't be used anyway
      }
    }
    if (applyDiff) {
      double diffk = 0.125 * (deltaAB + deltaCD - deltaAC - deltaBD);
      assert(isNumber(diffk));
      diff += diffk;
      sprime[k] += 1;
    }
    delta_Vsminus1_Wp[k] = deltaAB; // in our out the if???

    // update VsHasChildren for next call
    VsHasChildren[k] = Vsminus1HasChildren[k] || hasChildren[Ws->node_index][k];
  }
  ParallelContext::sumDouble(diff);
  double Ls = Lsminus1 + diff;
  if (Ls > bestLs) {
    bestLs = Ls;
    bestRegraftNode = Vs;
    bestS = s;
  }
  // recursive call
  if (s >= maxRadius) {
    return;
  }
  if (Vs->back->next) {
    getBestSPRRecMissing(s+1, sprime, W0s, Wp, Ws, Vs, delta_Vsminus1_Wp, Vs->back->next, 
        Ls, bestRegraftNode, bestLs, bestS, subBMEs, belongsToPruned, hasChildren, Vsminus1HasChildren, VsHasChildren);
    getBestSPRRecMissing(s+1, sprime, W0s, Wp, Ws, Vs, delta_Vsminus1_Wp, Vs->back->next->next, 
        Ls, bestRegraftNode, bestLs, bestS, subBMEs, belongsToPruned, hasChildren, Vsminus1HasChildren, VsHasChildren);
  }
}


void MiniBME::getBestSPR(PLLUnrootedTree &speciesTree,
      pll_unode_t *&bestPruneNode,
      pll_unode_t *&bestRegraftNode,
      double &bestDiff)
{
  bestDiff = 0.0;
  unsigned int bestS = 1;
  for (auto pruneNode: speciesTree.getPostOrderNodes()) {
    if (getBestSPRFromPrune(pruneNode, bestRegraftNode, bestDiff, bestS)) {
      bestPruneNode = pruneNode;
    }
  }
  Logger::info << "Best S=" << bestS << std::endl;
}


bool MiniBME::getBestSPRFromPrune(pll_unode_t *prunedNode,
      pll_unode_t *&bestRegraftNode,
      double &bestDiff,
      unsigned int &bestS)
{
  bool foundBetter = false;
  double oldBestDiff = bestDiff;
  auto Wp = prunedNode;
  if (!Wp->back->next) {
    return foundBetter;
  }
  std::vector<pll_unode_t *> V0s;
  V0s.push_back(Wp->back->next->back);
  V0s.push_back(Wp->back->next->next->back);
  for (auto V0: V0s) {
    auto V = getOtherNext(Wp->back, V0->back);
    if (!V->back->next) {
      continue;
    }
    std::vector<pll_unode_t *> V1s;
    V1s.push_back(V->back->next);
    V1s.push_back(V->back->next->next);
    for (auto V1: V1s) {
      unsigned int s = 1;
      pll_unode_t *W0 = V0;
      pll_unode_t *Wsminus1 = W0;
      if (_prune) {
        unsigned int K = _subBMEs[0][0].size();
        std::vector<pll_unode_t *> W0s(K, nullptr);
        std::vector<double> delta_V0_Wp; // not used at first iteration
        std::vector<unsigned int> sprime(K, 1);
        std::vector<bool> V0HasChildren = _hasChildren[V0->node_index];
        std::vector<bool> Vminus1HasChildren(V0HasChildren.size()); // won't be read at first iteration
        getBestSPRRecMissing(s, sprime, W0s, Wp,  Wsminus1, V,delta_V0_Wp,
            V1, 0.0, bestRegraftNode, bestDiff, bestS, _subBMEs, 
            _belongsToPruned, _hasChildren, Vminus1HasChildren, V0HasChildren);
      } else {
        double deltaAB = 0.0; // not used at first iteration
        getBestSPRRec(s, W0, Wp,  Wsminus1, V, deltaAB,
            V1, 0.0, bestRegraftNode, bestDiff, _subBMEs);
      }
      if (bestDiff > oldBestDiff) {
        foundBetter = true;
      }
    }
  }
  return foundBetter;
}

