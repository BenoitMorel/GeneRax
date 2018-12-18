#include "UndatedDLModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

using namespace std;

UndatedDLModel::UndatedDLModel() :
  O_R(1),
  allCLVInvalid(true),
  dupRate_(0),
  lossRate_(0)
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
  vector<ScaledValue> zeros(getPrunedSpeciesCount());
  uq = vector<vector<ScaledValue> >(maxId + 1,zeros);
  ll = vector<ScaledValue>(getPrunedSpeciesCount());
  repeatsId = vector<unsigned long>(maxId + 1, 0);
  setRates(dupRate_, lossRate_);
}

double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

void UndatedDLModel::setRates(double dupRate, 
  double lossRate,
  double transferRate) {
  geneRoot_ = 0;
  dupRate_ = dupRate;
  lossRate_ = lossRate;
  if (!getPrunedSpeciesCount()) {
    return;
  }
  PD = vector<double>(getPrunedSpeciesCount(), dupRate);
  PL = vector<double>(getPrunedSpeciesCount(), lossRate);
  PS = vector<double>(getPrunedSpeciesCount(), 1.0);


  for (auto speciesNode: getPrunedSpecies()) {
    int e = getPrunedSpeciesIndex(speciesNode);
    double sum = PD[e] + PL[e] + PS[e];
    PD[e] /= sum;
    PL[e] /= sum;
    PS[e] /= sum;
  } 
  uE = vector<double>(getPrunedSpeciesCount(), 0.0);
  for (auto speciesNode: getPrunedSpecies()) {
    int e = getPrunedSpeciesIndex(speciesNode);
    double a = PD[e];
    double b = -1.0;
    double c = PL[e];
    if (speciesNode->left) {
      c += PS[e] * getExtProba(speciesNode->left)  * getExtProba(speciesNode->right);
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    ASSERT_PROBA(proba)
    getExtProba(speciesNode) = proba;
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
  if (!geneNode->next) {
    repeatsId[geneNode->node_index] = geneToSpecies_[geneNode->node_index]; 
  } else {
    long left = repeatsId[geneNode->next->back->node_index];
    long right = repeatsId[geneNode->next->next->back->node_index];
    repeatsId[geneNode->node_index] = min(right, left) + 42 * max(right, left);
  }
  for (auto speciesNode: getPrunedSpecies()) {
    computeProbability(geneNode, speciesNode, uq[geneNode->node_index][getPrunedSpeciesIndex(speciesNode)]);
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
  int e = getPrunedSpeciesIndex(speciesNode);
  int f = 0;
  int g = 0;
  if (!isSpeciesLeaf) {
    f = getPrunedSpeciesIndex(speciesNode->left);
    g = getPrunedSpeciesIndex(speciesNode->right);
  }
  proba = ScaledValue();
  if (isSpeciesLeaf and isGeneLeaf and speciesNode->node_index == geneToSpecies_[gid]) {
    // present
    proba = ScaledValue(PS[e], 0);
  }
  if (not isGeneLeaf) {
    int gp_i = leftGeneNode->node_index;
    int gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      proba += ScaledValue::superMult1(uq[gp_i][f], uq[gpp_i][g],
          uq[gp_i][g], uq[gpp_i][f],
          PS[e]);
    }
    // D event
    ScaledValue temp = uq[gp_i][e];
    temp *= uq[gpp_i][e];
    temp *= 2.0 * PD[e];
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
          ll[f], uE[g],
          ll[g], uE[f],
          PS[e]);
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
  for (auto speciesNode: getPrunedSpecies()) {
    int e = getPrunedSpeciesIndex(speciesNode);
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
  getIdsPostOrder(*treeinfo, geneIds); 

  // main loop
  updateCLVs(*treeinfo);
  
  /*
  unordered_set<unsigned long> plop;
  for (auto r: repeatsId) {
    plop.insert(r);
  } 
  Logger::info << plop.size() << "/" << geneIds.size() << endl;
*/
  ScaledValue bestValue;
  auto bestRoot = computeLikelihoods(*treeinfo, bestValue);
  if (bestRoot) {
    geneRoot_ = bestRoot;
  }
  ScaledValue total;
  if (!Arguments::rootedGeneTree) {
    for (int e = 0; e < getPrunedSpeciesCount(); e++) {
      total += ll[e];
    }
  } else {
    total = bestValue;
  }
  double res = total.getLogValue(); //log(root_sum / survive / O_norm * (speciesNodesCount_));
  return res;
}


