#include "UndatedDLModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

using namespace std;

UndatedDLModel::UndatedDLModel() :
  O_R(1)
{
  Logger::info << "creating undated dl model" << endl;
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
    uE[speciesNode->node_index] = (-b - sqrt(b * b - 4 * a * c)) / (2.0 * a);
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
  int gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = geneNode->next->back;
    rightGeneNode = geneNode->next->next->back;
  }
  for (auto speciesNode: speciesNodes_) {
    bool isSpeciesLeaf = !speciesNode->left;
    int e = speciesNode->node_index;
    int f = 0;
    int g = 0;
    int gp_i = 0;
    int gpp_i = 0;
    if (!isSpeciesLeaf) {
      f = speciesNode->left->node_index;
      g = speciesNode->right->node_index;
    }
    double uq_sum = 0;
    double uq_Sevent1 = 0;
    double uq_Sevent2 = 0;
    double uq_Devent = 0;
    double uq_SLevent1 = 0;
    double uq_SLevent2 = 0;
    if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
      // present
      uq_sum += PS[e];
    }
    if (not isGeneLeaf) {
      gp_i = leftGeneNode->node_index;
      gpp_i = rightGeneNode->node_index;
      if (not isSpeciesLeaf) {
        // S event
        uq_Sevent1 += PS[e] * uq[gp_i][f] * uq[gpp_i][g];
        uq_Sevent2 += PS[e] * uq[gp_i][g] * uq[gpp_i][f];
      }
      // D event
      uq_Devent += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2);
    }
    if (not isSpeciesLeaf) {
      // SL event
      uq_SLevent1 += PS[e] * uq[gid][f] * uE[g];
      uq_SLevent2 += PS[e] * uq[gid][g] * uE[f];
    }
    double uq_internal_sum = uq_Sevent1 + uq_Sevent2 + uq_Devent + uq_SLevent1 + uq_SLevent2;
    uq_sum += uq_internal_sum;
      
    for (int event = 0; event < (int)LAST; ++event) {
      eventsCount[gid][e][event] = 0;
    }
    if (uq_Sevent1 != 0) {
      uq_Sevent1 /= uq_internal_sum;
      for (int event = 0; event < (int)LAST; ++event) {
        eventsCount[gid][e][event] += uq_Sevent1 * eventsCount[gp_i][f][event];
      }
      eventsCount[gid][e][(int)SPEC] += uq_Sevent1;
    }
    if (uq_Sevent2 != 0) {
      uq_Sevent2 /= uq_internal_sum;
      for (int event = 0; event < (int)LAST; ++event) {
        eventsCount[gid][e][event] += uq_Sevent2 * eventsCount[gpp_i][g][event];
      }
      eventsCount[gid][e][(int)SPEC] += uq_Sevent2;
    }
    // DL event
    uq[gid][e] = uq_sum / (1.0 - 2.0 * PD[e] * uE[e]);
  }
}


pll_unode_t * UndatedDLModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo)
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, geneIds);
  pll_unode_t *bestRoot = 0;
  for (auto speciesNode: speciesNodes_) {
    bool isSpeciesLeaf = !speciesNode->left;
    int e = speciesNode->node_index;
    int f = 0;
    int g = 0;
    if (!isSpeciesLeaf) {
      f = speciesNode->left->node_index;
      g = speciesNode->right->node_index;
    }
    double uq_e = 0;
    for (auto root: roots) {
      double uq_sum = 0;
      int gp_i = root->node_index;
      int gpp_i = root->back->node_index;
      if (not isSpeciesLeaf) {
        uq_sum += PS[e] * (uq[gp_i][f] * uq[gpp_i][g] + uq[gp_i][g] * uq[gpp_i][f]);
      }
      // D event
      uq_sum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2);
      if (Arguments::rootedGeneTree) {
        if (uq_sum > uq_e) {
          uq_e = max(uq_sum, uq_e);
          bestRoot = root;
        }
      } else {
        uq_e += uq_sum;
      }
    }

    if (not isSpeciesLeaf) {
      // SL event
      uq_e += PS[e] * (ll[f] * uE[g] + ll[g] * uE[f]); 
    }
    ll[e] = uq_e / (1.0 - 2.0 * PD[e] * uE[e]);
  }
  return bestRoot;
}




double UndatedDLModel::computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
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
  
  vector<int> zerosCounts(LAST, 0);
  vector<vector <int> > zerosEvents(speciesNodesCount_, zerosCounts);
  eventsCount = vector<  vector<vector <int> > >(maxId + 1, zerosEvents);
  // main loop
  updateCLVs(*treeinfo);
  
  auto bestRoot = computeLikelihoods(*treeinfo);
  if (bestRoot) {
    geneRoot_ = bestRoot;
  }
  for (int e = 0; e < speciesNodesCount_; e++) {
    double O_p = 1;
    if (e == (speciesNodesCount_ - 1)) 
      O_p = O_R;
    O_norm += O_p;
    root_sum += ll[e] * O_p;
    survive += (1 - uE[e]);
  }
  return root_sum / survive / O_norm * (speciesNodesCount_);
}


