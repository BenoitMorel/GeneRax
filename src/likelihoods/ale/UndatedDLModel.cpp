#include "UndatedDLModel.hpp"
#include <Arguments.hpp>

using namespace std;

UndatedDLModel::UndatedDLModel() :
  O_R(1),
  geneRoot(0)
{
}



void UndatedDLModel::fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) 
{
  if (node->left) {
    assert(node->right);
    fillNodesPostOrder(node->left, nodes);
    fillNodesPostOrder(node->right, nodes);
  }
  nodes.push_back(node);
}


void UndatedDLModel::setSpeciesTree(pll_rtree_t *speciesTree)
{
  speciesNodesCount = speciesTree->tip_count + speciesTree->inner_count;
  speciesNodes.clear();
  fillNodesPostOrder(speciesTree->root, speciesNodes);
  speciesNameToId.clear();
  for (auto node: speciesNodes) {
    if (!node->left) {
      speciesNameToId[node->label] = node->node_index;
    }
  }
}

void UndatedDLModel::setRates(double dupRate, 
  double lossRate,
  double transferRates) {
  geneRoot = 0;
  PD = vector<double>(speciesNodesCount, dupRate);
  PL = vector<double>(speciesNodesCount, lossRate);
  PS = vector<double>(speciesNodesCount, 1.0);
  for (auto speciesNode: speciesNodes) {
    int e = speciesNode->node_index;
    double sum = PD[e] + PL[e] + PS[e];
    PD[e] /= sum;
    PL[e] /= sum;
    PS[e] /= sum;
  } 
  uE = vector<double>(speciesNodesCount, 0.0);
  for (auto speciesNode: speciesNodes) {
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


void getIdsPostOrderRec(pll_unode_t *node, 
    vector<bool> &marked,
    vector<int> &nodeIds)
{
  if (marked[node->node_index]) {
    return;
  }
  if (node->next) {
    getIdsPostOrderRec(node->next->back, marked, nodeIds);
    getIdsPostOrderRec(node->next->next->back, marked, nodeIds);
  }
  nodeIds.push_back(node->node_index);
  marked[node->node_index] = true;
}

void UndatedDLModel::getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds) {
  int nodesNumber = tree.subnode_count;
  nodeIds.clear();
  vector<bool> marked(nodesNumber, false);
  if (Arguments::aleRooted && geneRoot) {
    getIdsPostOrderRec(geneRoot, marked, nodeIds);
    getIdsPostOrderRec(geneRoot->back, marked, nodeIds);
    return;
  } 
  
  for (int i = 0; i < nodesNumber; ++i) {
    getIdsPostOrderRec(tree.subnodes[i], marked, nodeIds);
  }
}

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
  for (auto speciesNode: speciesNodes) {
    bool isSpeciesLeaf = !speciesNode->left;
    int e = speciesNode->node_index;
    int f = 0;
    int g = 0;
    if (!isSpeciesLeaf) {
      f = speciesNode->left->node_index;
      g = speciesNode->right->node_index;
    }
    double uq_sum = 0;
    if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies[gid]) {
      // present
      uq_sum += PS[e];
    }
    if (not isGeneLeaf) {
      int gp_i = leftGeneNode->node_index;
      int gpp_i = rightGeneNode->node_index;
      if (not isSpeciesLeaf) {
        uq_sum += PS[e] * (uq[gp_i][f] * uq[gpp_i][g] + uq[gp_i][g] * uq[gpp_i][f]);
      }
      // D event
      uq_sum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2);
    }
    if (not isSpeciesLeaf) {
      // SL event
      uq_sum += PS[e] * (uq[gid][f] * uE[g] + uq[gid][g] * uE[f]);
    }
    uq[gid][e] = uq_sum / (1.0 - 2.0 * PD[e] * uE[e]);
  }
}

void UndatedDLModel::getRoots(pllmod_treeinfo_t &treeinfo, vector<pll_unode_t *> &roots)
{
  roots.clear();
  if (Arguments::aleRooted && geneRoot) {
    roots.push_back(geneRoot);
    return;
  }
  vector<bool> marked(geneIds.size(), false);
  for (auto id: geneIds) {
    auto node = treeinfo.subnodes[id];
    if (marked[node->clv_index] || marked[node->back->clv_index]) {
      continue;
    }
    roots.push_back(node);
    marked[node->clv_index] = true;
  }

}

pll_unode_t * UndatedDLModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo)
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots);
  pll_unode_t *bestRoot = 0;
  for (auto speciesNode: speciesNodes) {
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
      if (Arguments::aleRooted) {
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

void UndatedDLModel::mapGenesToSpecies(pllmod_treeinfo_t &treeinfo)
{
  geneToSpecies.resize(treeinfo.subnode_count);
  for (int i = 0; i < treeinfo.subnode_count; ++i) {
    auto node = treeinfo.subnodes[i];
    if (!node->next) {
      string speciesName = geneNameToSpeciesName[string(node->label)]; 
      geneToSpecies[node->node_index] = speciesNameToId[speciesName];
    }
  }
}


void UndatedDLModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  mapGenesToSpecies(*treeinfo);
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
  vector<double> zeros(speciesNodesCount, 0.0);
  uq = vector<vector<double>>(maxId + 1,zeros);
  ll = vector<double>(speciesNodesCount, 0.0);

  // main loop
  updateCLVs(*treeinfo);
  
  auto bestRoot = computeLikelihoods(*treeinfo);
  if (bestRoot) {
    geneRoot = bestRoot;
  }
  for (int e = 0; e < speciesNodesCount; e++) {
    double O_p = 1;
    if (e == (speciesNodesCount - 1)) 
      O_p = O_R;
    O_norm += O_p;
    root_sum += ll[e] * O_p;
    survive += (1 - uE[e]);
  }
  return root_sum / survive / O_norm * (speciesNodesCount);
}

void UndatedDLModel::setGeneSpeciesMap(const GeneSpeciesMapping &map)
{
  geneNameToSpeciesName = map.getMap();
}

