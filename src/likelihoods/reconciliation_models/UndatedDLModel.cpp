#include "UndatedDLModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

using namespace std;
const int CACHE_SIZE = 100000;

UndatedDLModel::UndatedDLModel():
  allCLVInvalid(true)
{
  maxId = 1;
  Logger::info << "creating undated dl model" << endl;
}

#define IS_PROBA(x) ((x) >= 0 && (x) <= 1 && !isnan(x))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));

void UndatedDLModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  AbstractReconciliationModel::setInitialGeneTree(treeinfo);
  getIdsPostOrder(*treeinfo, geneIds); 
  maxId = 0;
  for (auto gid: geneIds)
    maxId = max(maxId, gid);
  // init ua with zeros
  vector<ScaledValue> zeros(speciesNodesCount_);
  uq = vector<vector<ScaledValue> >(maxId + 1,zeros);
  virtual_uq = vector<vector<ScaledValue> > (maxId + 1,zeros);
  repeatsId = vector<unsigned long>(maxId + 1, 0);
  invalidateAllCLVs();
}

static double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

void UndatedDLModel::setRates(double dupRate, 
  double lossRate,
  double transferRates) {
  geneRoot_ = 0;
  cache_.clear();
  cache_.resize(CACHE_SIZE);
  PD = vector<double>(speciesNodesCount_, dupRate);
  PL = vector<double>(speciesNodesCount_, lossRate);
  PS = vector<double>(speciesNodesCount_, 1.0);
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    double sum = PD[e] + PL[e] + PS[e];
    PD[e] /= sum;
    PL[e] /= sum;
    PS[e] /= sum;
  } 
  uE = vector<double>(speciesNodesCount_, 0.0);
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    double a = PD[e];
    double b = -1.0;
    double c = PL[e];
    if (speciesNode->left) {
      c += PS[e] * uE[speciesNode->left->node_index]  * uE[speciesNode->right->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    ASSERT_PROBA(proba)
    uE[speciesNode->node_index] = proba;
  }
  invalidateAllCLVs();
}

UndatedDLModel::~UndatedDLModel() { }

void UndatedDLModel::markInvalidatedNodesRec(pll_unode_t *node)
{
  isCLVUpdated[node->node_index] = false;
  if (node->back->next) {
    markInvalidatedNodesRec(node->back->next);
    markInvalidatedNodesRec(node->back->next->next);
  }
}

void UndatedDLModel::markInvalidatedNodes(pllmod_treeinfo_t &treeinfo)
{
  for (int nodeIndex: invalidatedNodes) {
    auto node = treeinfo.subnodes[nodeIndex];
    markInvalidatedNodesRec(node);
  }
  invalidatedNodes.clear();
}

void UndatedDLModel::updateCLVsRec(pll_unode_t *node)
{
  if (isCLVUpdated[node->node_index]) {
    return;
  }
  if (node->next) {
    updateCLVsRec(node->next->back);
    updateCLVsRec(node->next->next->back);
  }
  updateCLV(node);
  isCLVUpdated[node->node_index] = true;
}

void UndatedDLModel::updateCLVs(pllmod_treeinfo_t &treeinfo)
{
  markInvalidatedNodes(treeinfo);
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, geneIds);
  for (auto root: roots) {
    updateCLVsRec(root);
    updateCLVsRec(root->back);
  }
}
void UndatedDLModel::computeGeneProbabilities(pll_unode_t *geneNode,
  vector<ScaledValue> &clv)
{
  for (auto speciesNode: speciesNodes_) {
    computeProbability(geneNode, 
        speciesNode, 
        clv[speciesNode->node_index]);
  }
}


void UndatedDLModel::updateCLV(pll_unode_t *geneNode)
{
#ifndef REPEATS
  computeGeneProbabilities(geneNode, uq[geneNode->node_index]);
#else
  int repeatId = 0;
  if (!geneNode->next) {
    repeatId = geneToSpecies_[geneNode->node_index] + 1;
  } else {
    long left = repeatsId[getLeft(geneNode, false)->node_index];
    long right = repeatsId[getRight(geneNode, false)->node_index];
    if (left * right) {
      repeatId = min(right, left) + (speciesNodesCount_+1) * max(right, left);// todobenoit this is wrong
    }
  }
  if (repeatId >= cache_.size()) {
    repeatId = 0;
  }
  repeatsId[geneNode->node_index] = repeatId;
  if (repeatId) {
    auto &cachedCLV = cache_[repeatId];
    if (!cachedCLV.size()) {
      computeGeneProbabilities(geneNode, uq[geneNode->node_index]);
      cachedCLV = uq[geneNode->node_index];
    } else {
      uq[geneNode->node_index] = cachedCLV;
    }
  } else {
    computeGeneProbabilities(geneNode, uq[geneNode->node_index]);
  }
#endif
}

void UndatedDLModel::invalidateCLV(int nodeIndex)
{
  invalidatedNodes.insert(nodeIndex);
}
  
void UndatedDLModel::invalidateAllCLVs()
{
  isCLVUpdated = vector<bool>(maxId + 1, false);
}

void UndatedDLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot) const
{
  int gid = isVirtualRoot ? geneNode->next->node_index : geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = getLeft(geneNode, isVirtualRoot);
    rightGeneNode = getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !speciesNode->left;
  int e = speciesNode->node_index;
  int f = 0;
  int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  proba = ScaledValue();
  if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
    // present
    proba = ScaledValue(PS[e], 0);
    return;
  }
  if (not isGeneLeaf) {
    int gp_i = leftGeneNode->node_index;
    int gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      proba += ScaledValue::superMult1(uq[gp_i][f], uq[gpp_i][g],
          uq[gp_i][g], uq[gpp_i][f],
          PS[e]);
    }
    // D event
    ScaledValue temp = uq[gp_i][e];
    temp *= uq[gpp_i][e];
    temp *= PD[e];
    proba += temp;
  }
  if (not isSpeciesLeaf) {
    // SL event
    if (!isVirtualRoot) {
      proba += ScaledValue::superMult2(
          uq[gid][f], uE[g],
          uq[gid][g], uE[f],
          PS[e]);
    } else {
      proba += ScaledValue::superMult2(
          virtual_uq[gid][f], uE[g],
          virtual_uq[gid][g], uE[f],
          PS[e]);
    }
  }
  // DL event
  proba /= (1.0 - 2.0 * PD[e] * uE[e]); 
  assert(proba.isProba());
}

void UndatedDLModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo)
{
  vector<ScaledValue> zeros(speciesNodesCount_); 
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, geneIds);
  for (auto root: roots) {
    int u = root->node_index;
    int uprime = root->back->node_index;
    for (auto speciesNode: speciesNodes_) {
      int e = speciesNode->node_index;
      pll_unode_t virtual_root;
      virtual_root.next = root;
      computeProbability(&virtual_root, speciesNode, virtual_uq[u][e], true);
      virtual_uq[uprime][e] = virtual_uq[u][e];
    }
  }
}

void UndatedDLModel::updateRoot(pllmod_treeinfo_t &treeinfo) 
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, geneIds);
  // find the best root
  ScaledValue max;
  for (auto root: roots) {
    ScaledValue sum;
    int u = root->node_index;
    for (auto speciesNode: speciesNodes_) {
      int e = speciesNode->node_index;
      sum += virtual_uq[u][e];
    }
    if (max < sum) {
      setRoot(root);
      max = sum;
    }
  }
}
double UndatedDLModel::getSumLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  ScaledValue total;
  vector<pll_unode_t *> roots;
  getRoots(*treeinfo, roots, geneIds);
  for (auto root: roots) {
    int u = root->node_index;
    ScaledValue sum;
    for (auto speciesNode: speciesNodes_) {
      int e = speciesNode->node_index;
      total += virtual_uq[u][e];
      sum += virtual_uq[u][e];
    }
  }
  double res = total.getLogValue(); //log(root_sum / survive / O_norm * (speciesNodesCount_));
  return res;
}

double UndatedDLModel::computeLogLikelihoodInternal(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  auto root = getRoot();
  getIdsPostOrder(*treeinfo, geneIds); 
  updateCLVs(*treeinfo);
  computeLikelihoods(*treeinfo);
  if (Arguments::rootedGeneTree) {// && !getRoot()) {
    updateRoot(*treeinfo);
    while (root != getRoot()) {
      updateCLVs(*treeinfo);
      computeLikelihoods(*treeinfo);
      root = getRoot();
      updateRoot(*treeinfo);
    }
  }
  return getSumLikelihood(treeinfo);
}


