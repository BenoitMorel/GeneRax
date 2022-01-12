#include "MiniBMEPruned.hpp"
#include <DistanceMethods/MiniNJ.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <search/UNNISearch.hpp>
#include <limits>

using DistanceVectorMatrix = std::vector<DistanceMatrix>;

static corax_unode_t *getOtherNext(corax_unode_t *n1, 
    corax_unode_t *n2)
{
  if (n1->next == n2) {
    return n2->next;
  } else {
    return n1->next;
  }
}

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


/**
 *  Computes the distances between a leaf (corresponding to
 *  currentNode at the first recursive call) and all other 
 *  nodes in the unrooted tree. If belongsToPruned is set,
 *  the distances are computed on the pruned subtree
 */
static void fillDistancesRec(corax_unode_t *currentNode, 
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

/**
 *  Fills the internode distance matrix of a pruned species tree.
 *  belongsToPruned tells which species leaves do not belong
 *  to the pruned species tree.
 */
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
    


MiniBMEPruned::MiniBMEPruned(const PLLUnrootedTree &speciesTree,
    const Families &families,
    double minbl)
{
  Logger::timed << "MiniBMEPruned::MiniBMEPruned 1" << std::endl;
  bool minMode = true;
  bool reweight = false;
  bool ustar = false;
  for (auto leaf: speciesTree.getLeaves()) {
    std::string label(leaf->label);
    _speciesStringToSpeciesId.insert({label, _speciesStringToSpeciesId.size()});
  }
  for (unsigned int i = 0; i < speciesTree.getLeavesNumber() + 3; ++i) {
    _pows.push_back(std::pow(0.5, i));
  }
  _N = speciesTree.getDirectedNodesNumber(); 
  // fill species-spid mappings
  Families perCoreFamilies;
  PerCoreGeneTrees::getPerCoreFamilies(families, perCoreFamilies);
  _K = perCoreFamilies.size();
  _geneDistanceMatrices.resize(_K);
  _geneDistanceDenominators.resize(_K);
  unsigned int speciesNumber = speciesTree.getLeavesNumber();
  _perFamilyCoverageStr.resize(_K);
  _perFamilyCoverage.resize(_K);
  // map each species string to a species ID
  std::set< std::vector<bool> > coverageSet;
  Logger::timed << "MiniBMEPruned::MiniBMEPruned 2" << std::endl;
  for (unsigned int i = 0; i < _K; ++i) {
    auto &family = perCoreFamilies[i];
    _perFamilyCoverage[i] = std::vector<bool>(speciesTree.getLeavesNumber(), false);
    _geneDistanceMatrices[i] = getNullMatrix(speciesNumber);
    _geneDistanceDenominators[i] = getNullMatrix(speciesNumber);
    GeneSpeciesMapping mappings;
    mappings.fill(family.mappingFile, family.startingGeneTree);
    for (auto species: mappings.getCoveredSpecies()) {
      _perFamilyCoverageStr[i].insert(species);
      _perFamilyCoverage[i][_speciesStringToSpeciesId.at(species)] = true;
    }
    coverageSet.insert(_perFamilyCoverage[i]);
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
          ustar,
          minbl);
    }
  }
  Logger::info << "Families: " << _K << std::endl;
  Logger::info << "Patterns: " << coverageSet.size() << std::endl;
  _prunedSpeciesMatrices = std::vector<DistanceMatrix>(
      _K, getNullMatrix(speciesNumber));
  auto nodesNumber = speciesTree.getDirectedNodesNumber();
  Logger::timed << "MiniBMEPruned::MiniBMEPruned 3" << std::endl;
  auto zero = getMatrix(nodesNumber, 
        _K, 
        std::numeric_limits<double>::infinity());
  Logger::timed << "MiniBMEPruned::MiniBMEPruned 4" << std::endl;
#ifdef OLDBME
  _subBMEs = DistanceVectorMatrix(nodesNumber);
  Logger::timed << "MiniBMEPruned::MiniBMEPruned 5" << std::endl;
  std::fill(_subBMEs.begin(), _subBMEs.end(), zero);
  Logger::timed << "MiniBMEPruned::MiniBMEPruned 6" << std::endl;
#else
  _subBMEs = std::vector<double>(_N * _N * _K, 
      std::numeric_limits<double>::infinity());
  // check integer overflow
  assert(_subBMEs.size() / _K == _N * _N);
#endif
  std::cerr << _subBMEs.size() << std::endl;
  for (auto sp1: speciesTree.getLeaves()) {
    auto i1 = sp1->node_index;
    for (auto sp2: speciesTree.getLeaves()) {
      auto i2 = sp2->node_index;
      for (auto k = 0; k < _K; ++k) {
        setCell(i1, i2, k, _geneDistanceMatrices[k][i1][i2]);
      }
    }
  }
  
  Logger::timed << "MiniBMEPruned::MiniBMEPruned 5" << std::endl;
}


double MiniBMEPruned::computeBME(const PLLUnrootedTree &speciesTree)
{
  unsigned int N = _speciesStringToSpeciesId.size();
  DistanceMatrix speciesDistanceMatrix = getNullMatrix(N);
  fillSpeciesDistances(speciesTree, 
      _speciesStringToSpeciesId,
      speciesDistanceMatrix);
  double res = 0.0;
  // O(kn^2)
  for (unsigned int k = 0; k < _K; ++k) {
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

        res += _geneDistanceMatrices[k][i][j] * _pows[_prunedSpeciesMatrices[k][i][j]];
      }
    }
  }
  ParallelContext::sumDouble(res);
  _computeSubBMEsPrune(speciesTree);
  return res;
}



 void MiniBMEPruned::_computeSubBMEsPruneRec(corax_unode_t *n1,
    corax_unode_t *n2,
    BoolMatrix &treated)
{
  auto i1 = n1->node_index;
  auto i2 = n2->node_index;
  unsigned int K = _K;
  if (treated[i1][i2]) {
    // stop if we've already been there
    return;
  }
  if (!n1->next && !n2->next) {
    // leaf-leaf case: already computed in constructor
  } else if (!n2->next) {
    // only n2 is a leaf, compute the symetric
    _computeSubBMEsPruneRec(n2, n1, treated); 
    for (unsigned int k = 0; k < _K; ++k) {
      setCell(i1, i2, k, getCell(i2, i1, k));
    }
  } else {
    // n2 is not a leaf, we can run recursion on n2
    auto left2 = n2->next->back;
    auto right2 = n2->next->next->back;
    _computeSubBMEsPruneRec(n1, left2, treated);
    _computeSubBMEsPruneRec(n1, right2, treated);
    for (unsigned int k = 0; k < K; ++k) {
      if (!_hasChildren[i1][k] || !_hasChildren[i2][k]) {
        // Edge case: either n1 or n2 is not in the path
        // of the induced tree, so its subBME should be
        // ignored
        setCell(i1, i2, k, 0.0);
        continue;
      }
      // Now we can assume that both n1 and n2 are in 
      // the path of the induced tree, but they might not
      // be binary nodes in the induced tree
      if (!_belongsToPruned[i2][k]) {
        // Case 1: n2 is not in the induced tree
        // Go down to n2's child that has children
        if (_hasChildren[left2->node_index][k]) {
          setCell(i1, i2, k, getCell(i1, left2->node_index, k));
        } else {
          setCell(i1, i2, k, getCell(i1, right2->node_index, k));
        }
      } else if (!_belongsToPruned[i1][k]) {
        // Case 2: n2 is in the induced tree, n1 is not
        // Note that n1 cannot be a leaf at this point
        // go down to n1's child that has children
        auto lefti1 = n1->next->back->node_index;
        auto righti1 = n1->next->next->back->node_index;
        if (_hasChildren[lefti1][k]) {
          setCell(i1, i2, k, getCell(lefti1, i2, k));
        } else {
          setCell(i1, i2, k, getCell(righti1, i2, k));
        }
      } else {
        // case 3: both n1 and n2 are in the induced tree
        double v = getCell(i1, left2->node_index, k);
        v += getCell(i1, right2->node_index, k);
        v *= 0.5;
        setCell(i1, i2, k, v);
      }
  
    }
  }
  treated[i1][i2] = true;
}

void MiniBMEPruned::_computeSubBMEsPrune(const PLLUnrootedTree &speciesTree)
{
  auto nodesNumber = speciesTree.getDirectedNodesNumber();
  // Fill hasChildren and belongsToPruned
  _hasChildren = getBoolMatrix(nodesNumber, 
      _K,
      false);
  _belongsToPruned = _hasChildren;
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto index = node->node_index;
    for (unsigned int k = 0; k < _K; ++k) {
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
  /*
  // this is a code optimization that does not work with
  // applying multiple SPR rounds simultaneously
  if (false && _toUpdate.before) {
    auto beforeNodes = speciesTree.getPostOrderNodesFrom(_toUpdate.before);
    auto afterNodes = speciesTree.getPostOrderNodesFrom(_toUpdate.after);
    auto prunedNodes = speciesTree.getPostOrderNodesFrom(_toUpdate.pruned);
    auto betweenNodes = speciesTree.getPostOrderNodesFrom(_toUpdate.between[0]);
    for (auto between: _toUpdate.between) {
      for (auto node: speciesTree.getPostOrderNodesFrom(between)) {
        betweenNodes.push_back(node);
      }
    }
    for (auto i1: prunedNodes) {
      for (auto i2: prunedNodes) {
        treated[i1->node_index][i2->node_index] = true;
      }
    }
    for (auto i1: betweenNodes) {
      for (auto i2: betweenNodes) {
        treated[i1->node_index][i2->node_index] = true;
      }
    }
    for (auto i1: afterNodes) {
      for (auto i2: afterNodes) {
        treated[i1->node_index][i2->node_index] = true;
      }
    }
    for (auto i1: beforeNodes) {
      for (auto i2: beforeNodes) {
        treated[i1->node_index][i2->node_index] = true;
      }
    }
    for (auto i1: beforeNodes) {
      for (auto i2: afterNodes) {
        treated[i1->node_index][i2->node_index] = true;
      }
    }
  }
  */
  for (auto n1: speciesTree.getPostOrderNodes()) {
    // we only need the subBMEs of the nodes of the
    // subtree rooted at n1->back
    for (auto n2: speciesTree.getPostOrderNodesFrom(n1->back)) {
      if (n2 == n1->back) {
        continue;
      }
      _computeSubBMEsPruneRec(n1,
        n2,
        treated);
    }
  }
}


void MiniBMEPruned::_getBestSPRRecMissing(unsigned int s,
    StopCriterion stopCriterion,
    std::vector<unsigned int> sprime, 
    std::vector<corax_unode_t *> W0s, 
    corax_unode_t *Wp, 
    corax_unode_t *Wsminus1, 
    corax_unode_t *Vsminus1, 
    std::vector<double> delta_Vsminus2_Wp, // previous deltaAB
    corax_unode_t *Vs, 
    double Lsminus1, // L_s-1
    corax_unode_t *&bestRegraftNode,
    SubBMEToUpdate &subBMEToUpdate,
    double &bestLs,
    unsigned int &bestS,
    const BoolMatrix &belongsToPruned,
    const BoolMatrix &hasChildren,
    std::vector<bool> Vsminus2HasChildren, // does Vsminus2 have children after the previous moves
    std::vector<bool> Vsminus1HasChildren) // does Vsminus1 have children after the previous moves
{
  if (s > stopCriterion.maxRadius || 
      stopCriterion.noImprovement > stopCriterion.maxRadiusWithoutImprovement) {
    return;
  }
  subBMEToUpdate.tempBetween.push_back(Wsminus1);
  corax_unode_t *Ws = getOtherNext(Vs, Vsminus1->back)->back; 
  unsigned int K = _K;
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
    double deltaCD = getCell(Ws->node_index, Vs->back->node_index, k);
    double deltaBD = getCell(Wp->node_index, Vs->back->node_index, k);
    
    // deltaAB and deltaAC are affected by the previous moves
    double deltaAB = 0.0;
    double deltaAC = 0.0;
    if (sprime[k] <= 1 && applyDiff) {
      deltaAB = getCell(W0s[k]->node_index, Wp->node_index, k);
      deltaAC = getCell(W0s[k]->node_index, Ws->node_index, k);
    } else if (W0s[k]) {
      deltaAC = _pows[sprime[k]] * 
        (getCell(W0s[k]->node_index, Ws->node_index, k) - 
         getCell(Wp->node_index, Ws->node_index, k));
      deltaAC += getCell(Vsminus1->node_index, Ws->node_index, k);
      if (hasChildren[Wsminus1->node_index][k] && Vsminus2HasChildren[k]) {
        deltaAB = 0.5 * (delta_Vsminus2_Wp[k] + 
          getCell(Wsminus1->node_index, Wp->node_index, k));
      } else if (hasChildren[Wsminus1->node_index][k] && !Vsminus2HasChildren[k]) {
        deltaAB = getCell(Wsminus1->node_index, Wp->node_index, k);
      } else if (!hasChildren[Wsminus1->node_index][k] && Vsminus2HasChildren[k]) {
        deltaAB = delta_Vsminus2_Wp[k];
      } else {
        deltaAB = std::numeric_limits<double>::infinity(); // won't be used anyway
      }
    }
    if (applyDiff) {
      double diffk = 0.125 * (deltaAB + deltaCD - deltaAC - deltaBD);
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
    subBMEToUpdate.between = subBMEToUpdate.tempBetween;
    subBMEToUpdate.after = Vs->back;
    subBMEToUpdate.pruned = Wp;
    stopCriterion.noImprovement = 0;
  } else {
    stopCriterion.noImprovement++;
  }
  if (Vs->back->next) {
    _getBestSPRRecMissing(s+1, stopCriterion, sprime, W0s, Wp, Ws, Vs, delta_Vsminus1_Wp, Vs->back->next, 
        Ls, bestRegraftNode, subBMEToUpdate, bestLs, bestS, belongsToPruned, hasChildren, Vsminus1HasChildren, VsHasChildren);
    _getBestSPRRecMissing(s+1, stopCriterion, sprime, W0s, Wp, Ws, Vs, delta_Vsminus1_Wp, Vs->back->next->next, 
        Ls, bestRegraftNode, subBMEToUpdate, bestLs, bestS, belongsToPruned, hasChildren, Vsminus1HasChildren, VsHasChildren);
  }
  subBMEToUpdate.tempBetween.pop_back();
}


void MiniBMEPruned::getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves)
{
  _toUpdate.reset();
  for (auto pruneNode: speciesTree.getPostOrderNodes()) {
    double bestDiff = 0.0;
    unsigned int bestS = 1;
    corax_unode_t *bestRegraftNode = nullptr;
    if (getBestSPRFromPrune(maxRadiusWithoutImprovement, pruneNode, bestRegraftNode, bestDiff, bestS)) {
      bestMoves.push_back(SPRMove(pruneNode, bestRegraftNode, bestDiff));
    }
  }
  std::sort(bestMoves.begin(), bestMoves.end());
}



bool MiniBMEPruned::getBestSPRFromPrune(unsigned int maxRadiusWithoutImprovement,
      corax_unode_t *prunedNode,
      corax_unode_t *&bestRegraftNode,
      double &bestDiff,
      unsigned int &bestS)
{
  bool foundBetter = false;
  double oldBestDiff = bestDiff;
  auto Wp = prunedNode;
  if (!Wp->back->next) {
    return foundBetter;
  }
  std::vector<corax_unode_t *> V0s;
  V0s.push_back(Wp->back->next->back);
  V0s.push_back(Wp->back->next->next->back);
  unsigned int K = _K;
  std::vector<unsigned int> sprime(K, 1);
  StopCriterion stopCriterion;
  stopCriterion.maxRadiusWithoutImprovement = maxRadiusWithoutImprovement;
  for (auto V0: V0s) {
    auto V = getOtherNext(Wp->back, V0->back);
    if (!V->back->next) {
      continue;
    }
    std::vector<corax_unode_t *> V1s;
    V1s.push_back(V->back->next);
    V1s.push_back(V->back->next->next);
    for (auto V1: V1s) {
      unsigned int s = 1;
      corax_unode_t *W0 = V0;
      corax_unode_t *Wsminus1 = W0;
      std::vector<corax_unode_t *> W0s(K, nullptr);
      std::vector<double> delta_V0_Wp; // not used at first iteration
      std::vector<bool> V0HasChildren = _hasChildren[V0->node_index];
      std::vector<bool> Vminus1HasChildren(V0HasChildren.size()); // won't be read at first iteration
      _getBestSPRRecMissing(s, stopCriterion, sprime, W0s, Wp,  Wsminus1, V,delta_V0_Wp,
          V1, 0.0, bestRegraftNode, _toUpdate, bestDiff, bestS, 
          _belongsToPruned, _hasChildren, Vminus1HasChildren, V0HasChildren);
      if (bestDiff > oldBestDiff) {
        foundBetter = true;
        _toUpdate.before = V0;
        oldBestDiff = bestDiff;
      }
    }
  }
  return foundBetter;
}

