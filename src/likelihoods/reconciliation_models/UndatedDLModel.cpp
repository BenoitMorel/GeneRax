#include "UndatedDLModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

using namespace std;

UndatedDLModel::UndatedDLModel() :
  O_R(1),
  allCLVInvalid(true)
{
  Logger::info << "creating undated dl model" << endl;
}

#define IS_PROBA(x) ((x) >= 0 && (x) <= 1 && !isnan(x))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));

  
void UndatedDLModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  AbstractReconciliationModel::setInitialGeneTree(treeinfo);
  getIdsPostOrder(*treeinfo, geneIds); 
  int maxId = 0;
  for (auto gid: geneIds)
    maxId = max(maxId, gid);
  // init ua with zeros
  vector<ScaledValue> zeros(speciesNodesCount_);
  uq = vector<vector<ScaledValue> >(maxId + 1,zeros);
  ll = vector<ScaledValue>(speciesNodesCount_);
  vector<int> zerosInt(speciesNodesCount_, 0);
}

double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

void UndatedDLModel::setRates(double dupRate, 
  double lossRate,
  double transferRates) {
  geneRoot_ = 0;
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
  allCLVInvalid = true;
}

UndatedDLModel::~UndatedDLModel() { }

bool isPresent(pll_unode_t *node, const unordered_set<int> &set) 
{
  assert(node);
  return set.find(node->node_index) != set.end();
}

bool getCLVsToUpdateRec(pll_unode_t *node, 
  const unordered_set<int> &invalidCLVs, 
  unordered_set<int> &nodesToUpdate, 
  unordered_set<int> &marked) 
{
  assert(node);
  if (isPresent(node, marked)) {
    return isPresent(node, nodesToUpdate);
  }
  marked.insert(node->node_index);
  bool needToUpdate = false; 
  if (isPresent(node, invalidCLVs)) {
    needToUpdate = true;
  }
  if (node->next) {
    needToUpdate |= getCLVsToUpdateRec(node->next->back, invalidCLVs, nodesToUpdate, marked);
    needToUpdate |= getCLVsToUpdateRec(node->next->next->back, invalidCLVs, nodesToUpdate, marked);
  }
  if (needToUpdate) {
    nodesToUpdate.insert(node->node_index);
  }  
  return needToUpdate;
}

void UndatedDLModel::getCLVsToUpdate(pllmod_treeinfo_t &treeinfo, unordered_set<int> &nodesToUpdate)
{
  nodesToUpdate.clear();
  if (allCLVInvalid) {
    for (int i = 0; i < (int) geneIds.size(); i++) {
      nodesToUpdate.insert(geneIds[i]);
    }
  } else {
    unordered_set<int> marked;
    if (!getRoot()) {
      for (int i = 0; i < (int) geneIds.size(); i++) {
        auto node = treeinfo.subnodes[geneIds[i]];    
        assert(node);
        getCLVsToUpdateRec(node, invalidCLVs, nodesToUpdate, marked);
      }
    } else {
      getCLVsToUpdateRec(getRoot(), invalidCLVs, nodesToUpdate, marked);
      getCLVsToUpdateRec(getRoot()->back, invalidCLVs, nodesToUpdate, marked);
    }
  }
}

void UndatedDLModel::updateCLVs(pllmod_treeinfo_t &treeinfo)
{
  unordered_set<int>  nodesToUpdate;
  if (!allCLVInvalid) {
    getCLVsToUpdate(treeinfo, nodesToUpdate);
  }
  for (int i = 0; i < (int) geneIds.size(); i++) {
    auto node = treeinfo.subnodes[geneIds[i]];
    if (allCLVInvalid || isPresent(node, nodesToUpdate)) {
      updateCLV(node);
    }
  }
  allCLVInvalid = false;
  invalidCLVs.clear();
}

void UndatedDLModel::updateCLV(pll_unode_t *geneNode)
{
  for (auto speciesNode: speciesNodes_) {
    computeProbability(geneNode, speciesNode, uq[geneNode->node_index][speciesNode->node_index]);
  }
}

void UndatedDLModel::invalidateCLV(int nodeIndex)
{
  invalidCLVs.insert(nodeIndex);
}

void UndatedDLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot)
{
  int gid = geneNode->node_index;
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
  }
  if (not isGeneLeaf) {
    int gp_i = leftGeneNode->node_index;
    int gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      proba += uq[gp_i][f] * uq[gpp_i][g] * PS[e];
      proba += uq[gp_i][g] * uq[gpp_i][f] * PS[e]; 
    }
    // D event
    proba += uq[gp_i][e] * uq[gpp_i][e] * 2.0 * PD[e];
  }
  if (not isSpeciesLeaf) {
    // SL event
    if (!isVirtualRoot) {
      proba += uq[gid][f] * uE[g] * PS[e];
      proba += uq[gid][g] * uE[f] * PS[e];
    } else {
      proba += ll[f] * uE[g] * PS[e];
      proba += ll[g] * uE[f] * PS[e];
    }
  }
  proba /= (1.0 - 2.0 * PD[e] * uE[e]); 
}

pll_unode_t * UndatedDLModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo,
    ScaledValue &bestValue)
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, geneIds);
  pll_unode_t *bestRoot = 0;
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    ScaledValue total;
    for (auto root: roots) {
      pll_unode_t virtual_root;
      virtual_root.next = root;
      ScaledValue value;
      computeProbability(&virtual_root, speciesNode, value, true);
      if (Arguments::rootedGeneTree) {
        if (bestValue < value) {
          bestValue = value;
          bestRoot = root;
        }
      } else {
        total += value;
      }
    }
    ll[e] = total;
  }
  return bestRoot;
}




double UndatedDLModel::computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  double survive = 0;
  double root_sum = 0;
  double O_norm = 0;

  getIdsPostOrder(*treeinfo, geneIds); 

  // main loop
  updateCLVs(*treeinfo);
  ScaledValue bestValue;
  auto bestRoot = computeLikelihoods(*treeinfo, bestValue);
  if (bestRoot) {
    geneRoot_ = bestRoot;
  }
  ScaledValue total;
  if (!Arguments::rootedGeneTree) {
    for (int e = 0; e < speciesNodesCount_; e++) {
      //O_norm += 1;
      total += ll[e];
      //root_sum += ll[e];
      //survive += (1 - uE[e]);
    }
  } else {
    total = bestValue;
  }
  double res = total.getLogValue(); //log(root_sum / survive / O_norm * (speciesNodesCount_));
  return res;
}


