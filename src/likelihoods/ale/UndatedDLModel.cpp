#include "UndatedDLModel.hpp"

// Bpp includes
#include <Bpp/BppString.h>
#include <Bpp/Numeric/Random/RandomTools.h>

#include <ale/tools/IO/IO.h>
#include <ale/tools/PhyloTreeToolBox.h>

using namespace bpp;
using namespace std;

UndatedDLModel::UndatedDLModel() :
  O_R(1)
{
}



void fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) 
{
  if (node->left) {
    assert(node->right);
    fillNodesPostOrder(node->left, nodes);
    fillNodesPostOrder(node->right, nodes);
  }
  nodes.push_back(node);
}


void UndatedDLModel::setSpeciesTree(const string &Sstring, pll_rtree_t *speciesTree)
{
  daughter.clear();
  son.clear();
  map<string, shared_ptr<bpp::PhyloNode>> name_node;
  map<shared_ptr<bpp::PhyloNode>, string> node_name;
  map<shared_ptr<bpp::PhyloNode>, int> node_ids;                         //Map between node and its id.
  S = shared_ptr<bpp::PhyloTree>(IO::newickToPhyloTree(Sstring, true));
  // For each node from the speciestree
  //vector<shared_ptr<bpp::PhyloNode>> nodes = S->getAllNodes();
  auto nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalRecursive_list(*S);

  // Set each branch length to 1.
  for (auto it = nodes.begin(); it != nodes.end(); it++) {
    if (S->hasFather(*it)) {
      auto edge_to_father = S->getEdgeToFather(*it);
      if (edge_to_father) S->getEdgeToFather(*it)->setLength(1);
    }
  }

  // For each node from the speciestree, record node names.
  // name_node gives a node according to the name
  // node_name gives a name accoring to the node.
  // * For leaves, get name.
  // * For internal nodes, gives a node name according to leaves under.
  // \todo node_name is useless until PhyloNode contains a name.
  for (auto it = nodes.begin(); it != nodes.end(); it++) {
    if (S->isLeaf(*it)) {
      name_node[(*it)->getName()] = (*it);
      node_name[(*it)] = (*it)->getName();
    } else {
      vector <string> leafnames = PhyloTreeToolBox::getLeavesNames(*S, (*it));
      sort(leafnames.begin(), leafnames.end());
      stringstream name;
      for (auto leafname = leafnames.begin(); leafname != leafnames.end(); leafname++)
        name << (*leafname) << ".";

      name_node[name.str()] = (*it);
      node_name[(*it)] = name.str();
    }
  }
  daughter.resize(nodes.size());
  son.resize(nodes.size());

  // register species
  last_branch = 0;
  speciesLastLeaf = 0;

  // associate each node to an id (node_ids and i_nodes).
  // this will determines the number of leaves (speciesLastLeaf)
  // and starting to determine the number of branches.
  // \todo Use PhyloTreeToolBox to get direct post-order.
  set<shared_ptr<bpp::PhyloNode>> saw;
  
  for (auto it = name_node.begin(); it != name_node.end(); it++)
    if (S->isLeaf((*it).second)) {
      auto node = (*it).second;
      node_ids[node] = last_branch;
      last_branch++;
      speciesLastLeaf++;
      saw.insert(node);
      // a leaf
      daughter[last_branch] = -1;
      // a leaf
      son[last_branch] = -1;
    }

  //ad-hoc postorder, each internal node is associated to an id and an id to an internal node
  //During the node exploration, set "ID" property to each node with value as the name.
  //
  vector<shared_ptr<bpp::PhyloNode>> next_generation;
  for (auto it = name_node.begin(); it != name_node.end(); it++)
    if (S->isLeaf((*it).second)) {
      auto node = (*it).second;
      next_generation.push_back(node);
    }

  while (next_generation.size()) {
    vector<shared_ptr<bpp::PhyloNode>> new_generation;
    for (auto it = next_generation.begin(); it != next_generation.end(); it++) {
      auto node = (*it);
      if (S->hasFather(node)) {
        auto father = S->getFather(node);
        auto sons = S->getSons(father);
        decltype(node) sister;

        if (sons[0] == node) 
          sister = sons[1];
        else 
          sister = sons[0];

        if (not node_ids.count(father) and saw.count(sister)) {
          node_ids[father] = last_branch;
          last_branch++;
          saw.insert(father);
          new_generation.push_back(father);
        }
      }
    }
    next_generation.clear();
    for (auto it = new_generation.begin();
         it != new_generation.end(); it++)
      next_generation.push_back((*it));
  }

  // Fill daughter and son maps. Daughter contains son left and son son right of a node.
  for (auto it = name_node.begin(); it != name_node.end(); it++) {
    if (not S->isLeaf(it->second))//not (*it).second->isLeaf())
    {
      auto node = (*it).second;
      auto sons = S->getSons(node);//node->getSons();
      daughter[node_ids[node]] = node_ids[sons[0]];
      son[node_ids[node]] = node_ids[sons[1]];
    }
  }
  
  speciesNodes.clear();
  fillNodesPostOrder(speciesTree->root, speciesNodes);
  speciesNameToLibpllId.clear();
  for (auto node: speciesNodes) {
    if (!node->left) {
      speciesNameToLibpllId[node->label] = node->node_index;
    }
  }

  
  for (auto node: nodes) {
    if (S->isLeaf(node)) {
      speciesNameToId[node->getName()] = node_ids[node];
    }
  }
}

void UndatedDLModel::setRates(double dupRate, double lossRate) {
  PD = dupRate;
  PL = lossRate;
  PS = 1;
  double sum = PD + PL + PS;
  PD /= sum;
  PL /= sum;
  PS /= sum;
  uE = vector<double>(last_branch, 0.0);
  for (int e = 0; e < speciesLastLeaf; ++e) {
    double a = PD;
    double b = -1.0;
    double c = PL;
    uE[e] = (-b - sqrt(b * b - 4 * a * c)) / (2.0 * a);
  }
  for (int e = speciesLastLeaf; e < last_branch; ++e) {
    int f = daughter[e];
    int g = son[e];
    double a = PD;
    double b = -1.0;
    double c = PL + PS * uE[f]  * uE[g];
    uE[e] = (-b - sqrt(b * b - 4 * a * c)) / (2.0 * a);
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

void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds) {
  int nodesNumber = tree.subnode_count;
  nodeIds.clear();
  vector<bool> marked(nodesNumber, false);
  for (int i = 0; i < nodesNumber; ++i) {
    getIdsPostOrderRec(tree.subnodes[i], marked, nodeIds);
  }
}

void UndatedDLModel::updateCLVs(pllmod_treeinfo_t &treeinfo)
{
  for (int i = 0; i < (int) geneIds.size(); i++) {
    auto gid = geneIds[i];
    pll_unode_t *geneNode = treeinfo.subnodes[gid];
    pll_unode_t *leftGeneNode = 0;     
    pll_unode_t *rightGeneNode = 0;     
    bool isGeneLeaf = !geneNode->next;
    if (!isGeneLeaf) {
      leftGeneNode = geneNode->next->back;
      rightGeneNode = geneNode->next->next->back;
    }
    for (int e = 0; e < last_branch; e++) {
      bool s_is_leaf = (e < speciesLastLeaf);
      int f = 0;
      int g = 0;
      if (not s_is_leaf) {
        f = daughter[e];
        g = son[e];
      }
      double uq_sum = 0;
      if (e < speciesLastLeaf and isGeneLeaf and e == geneToSpecies[gid]) {
        // present
        uq_sum += PS;
      }
      if (not isGeneLeaf) {
        int gp_i = leftGeneNode->node_index;
        int gpp_i = rightGeneNode->node_index;
        if (not s_is_leaf) {
          uq_sum += PS * (uq[gp_i][f] * uq[gpp_i][g] + uq[gp_i][g] * uq[gpp_i][f]);
        }
        // D event
        uq_sum += PD * (uq[gp_i][e] * uq[gpp_i][e] * 2);
      }
      if (not s_is_leaf) {
        // SL event
        uq_sum += PS * (uq[gid][f] * uE[g] + uq[gid][g] * uE[f]);
      }
      uq[gid][e] = uq_sum / (1.0 - 2.0 * PD * uE[e]);
    }
  }
}

void UndatedDLModel::getRoots(pllmod_treeinfo_t &treeinfo, vector<pll_unode_t *> &roots)
{
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

void UndatedDLModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo)
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots);

  for (int e = 0; e < last_branch; e++) {
    bool s_is_leaf = (e < speciesLastLeaf);
    int f = 0;
    int g = 0;
    if (not s_is_leaf) {
      f = daughter[e];
      g = son[e];
    }
    double uq_sum = 0;
    for (auto root: roots) {
      int gp_i = root->node_index;
      int gpp_i = root->back->node_index;
      if (not s_is_leaf) {
        uq_sum += PS * (uq[gp_i][f] * uq[gpp_i][g] + uq[gp_i][g] * uq[gpp_i][f]);
      }
      // D event
      uq_sum += PD * (uq[gp_i][e] * uq[gpp_i][e] * 2);
    }
    if (not s_is_leaf) {
      // SL event
      uq_sum += PS * (ll[f] * uE[g] + ll[g] * uE[f]); //todobenoit I am not sure why ll here
    }
    ll[e] = uq_sum / (1.0 - 2.0 * PD * uE[e]);
  }
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


double UndatedDLModel::pun(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  double survive = 0;
  double root_sum = 0;
  double O_norm = 0;

  // init gene ids
  
  getIdsPostOrder(*treeinfo, geneIds);
  inverseGeneIds.resize(geneIds.size());
  mapGenesToSpecies(*treeinfo);
  for (int i = 0; i < geneIds.size(); ++i) {
    inverseGeneIds[geneIds[i]] = i;
  }
 
  // init ua with zeros
  vector<double> zeros(last_branch, 0.0);
  uq = vector<vector<double>>(geneIds.size(),zeros);
  ll = vector<double>(last_branch, 0.0);

  // main loop
  updateCLVs(*treeinfo);
  computeLikelihoods(*treeinfo);
  for (int e = 0; e < last_branch; e++) {
    double O_p = 1;
    if (e == (last_branch - 1)) 
      O_p = O_R;
    O_norm += O_p;
    root_sum += ll[e] * O_p;
    survive += (1 - uE[e]);
  }
  return root_sum / survive / O_norm * (last_branch);
}

void UndatedDLModel::setMap(const GeneMap<string, string> &geneMap) { 
  geneNameToSpeciesName.clear();
  auto species = geneMap.getSpecies();
  for (auto s: species) {
    auto genes = geneMap.getGenes(s);
    for (auto g: genes) {
      geneNameToSpeciesName[g] = s; 
    }
  }
}

