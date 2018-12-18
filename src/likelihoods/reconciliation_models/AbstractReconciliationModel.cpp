#include "AbstractReconciliationModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

AbstractReconciliationModel::AbstractReconciliationModel():
  geneRoot_(0)
{

}

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

void AbstractReconciliationModel::getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds) {
  int nodesNumber = tree.subnode_count;
  nodeIds.clear();
  vector<bool> marked(nodesNumber, false);
  if (Arguments::rootedGeneTree && geneRoot_) {
    getIdsPostOrderRec(geneRoot_, marked, nodeIds);
    getIdsPostOrderRec(geneRoot_->back, marked, nodeIds);
    return;
  } 
  
  for (int i = 0; i < nodesNumber; ++i) {
    getIdsPostOrderRec(tree.subnodes[i], marked, nodeIds);
  }
}

void AbstractReconciliationModel::setGeneSpeciesMap(const GeneSpeciesMapping &map)
{
  geneNameToSpeciesName_ = map.getMap();
}

void AbstractReconciliationModel::mapGenesToSpecies(pllmod_treeinfo_t &treeinfo)
{
  vector<bool> speciesPresent(speciesNodesCount_, false);
  vector<bool> speciesToKeep(speciesNodesCount_, false);
  geneToSpecies_.resize(treeinfo.subnode_count);
  speciesToPrunedSpecies_.resize(speciesNodesCount_);
  for (int i = 0; i < treeinfo.subnode_count; ++i) {
    auto node = treeinfo.subnodes[i];
    if (!node->next) {
      string speciesName = geneNameToSpeciesName_[string(node->label)]; 
      geneToSpecies_[node->node_index] = speciesNameToId_[speciesName];
      speciesPresent[geneToSpecies_[node->node_index]] = true;
      speciesToKeep[geneToSpecies_[node->node_index]] = true;
    }
  }
  for (int e = 0; e < speciesNodesCount_; e++) {
    auto node = speciesNodes_[e];
    if (node->left) {
      speciesPresent[node->node_index] = speciesPresent[node->left->node_index] || speciesPresent[node->right->node_index];
      speciesToKeep[node->node_index] = speciesPresent[node->left->node_index] && speciesPresent[node->right->node_index];
    }

    // prunedSpeciesNodes_
    if ( speciesToKeep[node->node_index]) {
      speciesToPrunedSpecies_[node->node_index] = prunedSpeciesNodes_.size();
      prunedSpeciesNodes_.push_back(node);
    }
  }
  
  for (int i = 0; i < getPrunedSpeciesCount(); ++i) {
    assert(getPrunedSpeciesIndex(getPrunedSpecies()[i]) == i);
  }
  vector<bool> marked(getPrunedSpeciesCount(), false);
  for (auto speciesNode: getPrunedSpecies()) {
    if (speciesNode->left) {
      assert(marked[getPrunedSpeciesIndex(speciesNode->left)]);
      assert(marked[getPrunedSpeciesIndex(speciesNode->right)]);
    }
    assert(!marked[getPrunedSpeciesIndex(speciesNode)]);
    marked[getPrunedSpeciesIndex(speciesNode)] = true;
  }
  Logger::info << "pruned species nodes " << prunedSpeciesNodes_.size() << endl;
  Logger::info << "species nodes " << speciesNodesCount_ << endl;
}

void AbstractReconciliationModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  mapGenesToSpecies(*treeinfo);
}

void AbstractReconciliationModel::fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) 
{
  if (node->left) {
    assert(node->right);
    fillNodesPostOrder(node->left, nodes);
    fillNodesPostOrder(node->right, nodes);
  }
  nodes.push_back(node);
}


void AbstractReconciliationModel::setSpeciesTree(pll_rtree_t *speciesTree)
{
  speciesTree_ = speciesTree;
  speciesNodesCount_ = speciesTree->tip_count + speciesTree->inner_count;
  speciesNodes_.clear();
  fillNodesPostOrder(speciesTree->root, speciesNodes_);
  speciesNameToId_.clear();
  for (auto node: speciesNodes_) {
    if (!node->left) {
      speciesNameToId_[node->label] = node->node_index;
    }
  }
}

void AbstractReconciliationModel::getRoots(pllmod_treeinfo_t &treeinfo, 
    vector<pll_unode_t *> &roots,
    const vector<int> &geneIds)
{
  roots.clear();
  if (Arguments::rootedGeneTree && geneRoot_) {
    roots.push_back(geneRoot_);
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

pll_unode_t *AbstractReconciliationModel::getLeft(pll_unode_t *node, bool virtualRoot)
{
  return virtualRoot ? node->next : node->next->back;
}

pll_unode_t *AbstractReconciliationModel::getRight(pll_unode_t *node, bool virtualRoot)
{
  return virtualRoot ? node->next->back : node->next->next->back;
}

  
