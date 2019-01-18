#include "AbstractReconciliationModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

AbstractReconciliationModel::AbstractReconciliationModel():
  geneRoot_(0),
  firstCall_(true),
  _maxGeneId(1)
{
}
void AbstractReconciliationModel::init(pll_rtree_t *speciesTree, const GeneSpeciesMapping &map)
{
  setSpeciesTree(speciesTree);
  setGeneSpeciesMap(map);
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
  /*
  if (Arguments::rootedGeneTree && geneRoot_) {
    getIdsPostOrderRec(geneRoot_, marked, nodeIds);
    getIdsPostOrderRec(geneRoot_->back, marked, nodeIds);
    return;
  } 
  */
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
  geneToSpecies_.resize(treeinfo.subnode_count);
  for (int i = 0; i < treeinfo.subnode_count; ++i) {
    auto node = treeinfo.subnodes[i];
    if (!node->next) {
      string speciesName = geneNameToSpeciesName_[string(node->label)]; 
      geneToSpecies_[node->node_index] = speciesNameToId_[speciesName];
    }
  }
}

void AbstractReconciliationModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  mapGenesToSpecies(*treeinfo);
  getIdsPostOrder(*treeinfo, _geneIds); 
  _maxGeneId = 0;
  for (auto gid: _geneIds)
    _maxGeneId = max(_maxGeneId, gid);
  invalidateAllCLVs();
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
    if (geneRoot_->next) {
      roots.push_back(geneRoot_->next);
      roots.push_back(geneRoot_->next->next);
    }
    if (geneRoot_->back->next) {
      roots.push_back(geneRoot_->back->next);
      roots.push_back(geneRoot_->back->next->next);
    }
    return;
  }
  vector<bool> marked(geneIds.size(), false);
  for (auto id: geneIds) {
    auto node = treeinfo.subnodes[id];
    if (marked[node->node_index] || marked[node->back->node_index]) {
      continue;
    }
    roots.push_back(node->back);
    marked[node->node_index] = true;
  }
}
  
double AbstractReconciliationModel::computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  if (firstCall_) {
    setInitialGeneTree(treeinfo);
    firstCall_ = false;
  }
  return computeLogLikelihoodInternal(treeinfo);
}

pll_unode_t *AbstractReconciliationModel::getLeft(pll_unode_t *node, bool virtualRoot)
{
  return virtualRoot ? node->next : node->next->back;
}

pll_unode_t *AbstractReconciliationModel::getRight(pll_unode_t *node, bool virtualRoot)
{
  return virtualRoot ? node->next->back : node->next->next->back;
}

void AbstractReconciliationModel::markInvalidatedNodesRec(pll_unode_t *node)
{
  _isCLVUpdated[node->node_index] = false;
  if (node->back->next) {
    markInvalidatedNodesRec(node->back->next);
    markInvalidatedNodesRec(node->back->next->next);
  }
}

void AbstractReconciliationModel::markInvalidatedNodes(pllmod_treeinfo_t &treeinfo)
{
  for (int nodeIndex: _invalidatedNodes) {
    auto node = treeinfo.subnodes[nodeIndex];
    markInvalidatedNodesRec(node);
  }
  _invalidatedNodes.clear();
}

void AbstractReconciliationModel::updateCLVsRec(pll_unode_t *node)
{
  if (_isCLVUpdated[node->node_index]) {
    return;
  }
  if (node->next) {
    updateCLVsRec(node->next->back);
    updateCLVsRec(node->next->next->back);
  }
  updateCLV(node);
  _isCLVUpdated[node->node_index] = true;
}

void AbstractReconciliationModel::updateCLVs(pllmod_treeinfo_t &treeinfo)
{
  markInvalidatedNodes(treeinfo);
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, _geneIds);
  for (auto root: roots) {
    updateCLVsRec(root);
    updateCLVsRec(root->back);
  }
}

void AbstractReconciliationModel::invalidateCLV(int nodeIndex)
{
  _invalidatedNodes.insert(nodeIndex);
}
  
void AbstractReconciliationModel::invalidateAllCLVs()
{
  _isCLVUpdated = vector<bool>(_maxGeneId + 1, false);
}
  
