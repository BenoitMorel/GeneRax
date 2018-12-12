#include "UndatedDLModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

using namespace std;

UndatedDLModel::UndatedDLModel() :
  O_R(1)
{
  Logger::info << "creating undated dl model" << endl;
}

#define IS_PROBA(x) ((x) >= 0 && (x) <= 1 && !isnan(x))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));


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


}

UndatedDLModel::~UndatedDLModel() { }



void UndatedDLModel::updateCLVs(pllmod_treeinfo_t &treeinfo)
{
  for (int i = 0; i < (int) geneIds.size(); i++) {
    updateCLV(treeinfo.subnodes[geneIds[i]]);
  }
}

void UndatedDLModel::updateCLV(pll_unode_t *geneNode)
{
  for (auto speciesNode: speciesNodes_) {
    int scaler = 0;
    ScaledValue proba;
    computeProbability(geneNode, speciesNode, proba);
    uq[geneNode->node_index][speciesNode->node_index] = proba.value;
    uq_scalers[geneNode->node_index][speciesNode->node_index] = proba.scaler;
  }
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
  ListScaledValue total;
  if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
    // present
    total.add(PS[e], 0);
  }
  if (not isGeneLeaf) {
    int gp_i = leftGeneNode->node_index;
    int gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      total.add(PS[e] * uq[gp_i][f] * uq[gpp_i][g], 
          uq_scalers[gp_i][f] + uq_scalers[gpp_i][g]);
      total.add(PS[e] * uq[gp_i][g] * uq[gpp_i][f], 
          uq_scalers[gp_i][g] + uq_scalers[gpp_i][f]);
    }
    // D event
    total.add(PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2),
        uq_scalers[gp_i][e] + uq_scalers[gpp_i][e]);
  }
  if (not isSpeciesLeaf) {
    // SL event
    if (!isVirtualRoot) {
      total.add(PS[e] * uq[gid][f] * uE[g], uq_scalers[gid][f]);
      total.add(PS[e] * uq[gid][g] * uE[f], uq_scalers[gid][g]);
      //uq_sum += PS[e] * (uq[gid][f] * uE[g] + uq[gid][g] * uE[f]);
    } else {
      total.add(PS[e] * ll[f] * uE[g], ll_scalers[f]);
      total.add(PS[e] * ll[g] * uE[f], ll_scalers[g]);
      
      //uq_sum += PS[e] * (ll[f] * uE[g] + ll[g] * uE[f]);
    }
  }
  proba = total.total;
  proba.value /= (1.0 - 2.0 * PD[e] * uE[e]); 
}

pll_unode_t * UndatedDLModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo)
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, geneIds);
  pll_unode_t *bestRoot = 0;
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    ListScaledValue total;
    for (auto root: roots) {
      pll_unode_t virtual_root;
      virtual_root.next = root;
      ScaledValue value;
      computeProbability(&virtual_root, speciesNode, value, true);
      if (Arguments::rootedGeneTree) {
        assert(false);
        /*
        if (p > ll[e]) {
          ll[e] = p;
          bestRoot = root;
        }
        */
      } else {
        total.add(value);
      }
    }
    ll[e] = total.total.value;
    ll_scalers[e] = total.total.scaler;
  }
  return bestRoot;
}




double UndatedDLModel::computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  double survive = 0;
  double root_sum = 0;
  double O_norm = 0;

  getIdsPostOrder(*treeinfo, geneIds); 
  int maxId = 0;
  for (auto gid: geneIds)
    maxId = max(maxId, gid);
  // init ua with zeros
  vector<double> zeros(speciesNodesCount_, 0.0);
  uq = vector<vector<double>>(maxId + 1,zeros);
  ll = vector<double>(speciesNodesCount_, 0.0);
  
  vector<int> zerosInt(speciesNodesCount_, 0);
  uq_scalers = vector<vector<int>>(maxId + 1,zerosInt);
  ll_scalers = vector<int>(speciesNodesCount_, 0);

  // main loop
  updateCLVs(*treeinfo);
  
  auto bestRoot = computeLikelihoods(*treeinfo);
  if (bestRoot) {
    geneRoot_ = bestRoot;
  }
  ListScaledValue total;
  for (int e = 0; e < speciesNodesCount_; e++) {
    O_norm += 1;
    total.add(ll[e], ll_scalers[e]);
    //root_sum += ll[e];
    survive += (1 - uE[e]);
  }
  double res = total.total.getLogValue(); //log(root_sum / survive / O_norm * (speciesNodesCount_));
  return res;
}


