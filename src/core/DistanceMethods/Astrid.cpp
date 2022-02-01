#include "Astrid.hpp"
#include <DistanceMethods/MiniNJ.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
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
    




Astrid::Astrid(const PLLUnrootedTree &speciesTree,
    const Families &families,
    double minbl)
{
  bool minMode = true;
  bool reweight = false;
  bool ustar = false;
  for (auto leaf: speciesTree.getLeaves()) {
    std::string label(leaf->label);
    _speciesStringToSpeciesId.insert({label, _speciesIdToSpeciesString.size()});
    _speciesIdToSpeciesString.push_back(label);
  }
  for (unsigned int i = 0; i < speciesTree.getLeavesNumber() + 3; ++i) {
    _pows.push_back(std::pow(0.5, i));
  }
  _geneDistanceMatrices.resize(1);
  MiniNJ::computeDistanceMatrix(families,
    minMode, 
    reweight,
    ustar,
    minbl,
    _geneDistanceMatrices[0],
    _speciesIdToSpeciesString,
    _speciesStringToSpeciesId);
}


double Astrid::computeBME(const PLLUnrootedTree &speciesTree)
{
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
  _computeSubBMEs(speciesTree);
  return res;
}

static bool isNumber(double v) {
  return v == 0.0 || std::isnormal(v);
}


void Astrid::_computeSubBMEs(const PLLUnrootedTree &speciesTree)
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



static corax_unode_t *getOtherNext(corax_unode_t *n1, 
    corax_unode_t *n2)
{
  if (n1->next == n2) {
    return n2->next;
  } else {
    return n1->next;
  }
}

// test regrafting  Wp between Vs and Vsplus1 after 
// having regrafted Wp between Vsminus and Vs
// and recursively test further SPR mvoes
static void getBestSPRRec(unsigned int s,
    corax_unode_t *W0, 
    corax_unode_t *Wp, 
    corax_unode_t *Wsminus1, 
    corax_unode_t *Vsminus1, 
    double delta_Vsminus2_Wp, // previous deltaAB
    corax_unode_t *Vs, 
    double Lsminus1, // L_s-1
    corax_unode_t *&bestRegraftNode,
    double &bestLs,
    const std::vector<DistanceMatrix> &subBMEs)
{
  corax_unode_t *Ws = getOtherNext(Vs, Vsminus1->back)->back; 
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



void Astrid::getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves)
{
  for (auto pruneNode: speciesTree.getPostOrderNodes()) {
    double bestDiff = 0.0;
    unsigned int bestS = 1;
    corax_unode_t *bestRegraftNode = nullptr;
    if (getBestSPRFromPrune(pruneNode, bestRegraftNode, bestDiff, bestS)) {
      bestMoves.push_back(SPRMove(pruneNode, bestRegraftNode, bestDiff));
    }
  }
  std::sort(bestMoves.begin(), bestMoves.end());
}


bool Astrid::getBestSPRFromPrune(corax_unode_t *prunedNode,
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
      double deltaAB = 0.0; // not used at first iteration
      getBestSPRRec(s, W0, Wp,  Wsminus1, V, deltaAB,
          V1, 0.0, bestRegraftNode, bestDiff, _subBMEs);
      if (bestDiff > oldBestDiff) {
        foundBetter = true;
      }
    }
  }
  return foundBetter;
}

