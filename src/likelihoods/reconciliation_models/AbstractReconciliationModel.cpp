#include "AbstractReconciliationModel.hpp"
#include <IO/Logger.hpp>
  
AbstractReconciliationModel::AbstractReconciliationModel():
  rootedGeneTree_(false),
  geneRoot_(0),
  firstCall_(true),
  _maxGeneId(1)
{
}
void AbstractReconciliationModel::init(pll_rtree_t *speciesTree, const GeneSpeciesMapping &map, bool rootedGeneTree)
{
  rootedGeneTree_ = rootedGeneTree;
  setSpeciesTree(speciesTree);
  geneNameToSpeciesName_ = map.getMap();
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
  for (int i = 0; i < nodesNumber; ++i) {
    getIdsPostOrderRec(tree.subnodes[i], marked, nodeIds);
  }
}


void AbstractReconciliationModel::mapGenesToSpecies(pllmod_treeinfo_t &treeinfo)
{
  geneToSpecies_.resize(treeinfo.subnode_count);
  for (unsigned int i = 0; i < treeinfo.subnode_count; ++i) {
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
  if (rootedGeneTree_ && geneRoot_) {
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
  auto root = getRoot();
  updateCLVs(*treeinfo);
  computeLikelihoods(*treeinfo);
  if (rootedGeneTree_) {
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

void AbstractReconciliationModel::computeMLRoot(pllmod_treeinfo_t &treeinfo, 
    pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot) 
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, _geneIds);
  ScaledValue max;
  for (auto root: roots) {
    for (auto speciesNode: speciesNodes_) {
      ScaledValue ll = getRootLikelihood(treeinfo, root, speciesNode);
      if (max < ll) {
        max = ll;
        bestGeneRoot = root;
        bestSpeciesRoot = speciesNode;
      }
    }
  }
}

pll_unode_t *AbstractReconciliationModel::computeMLRoot(pllmod_treeinfo_t &treeinfo)
{
  pll_unode_t *bestRoot = 0;
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, _geneIds);
  ScaledValue max;
  for (auto root: roots) {
    ScaledValue rootProba = getRootLikelihood(treeinfo, root);
    if (max < rootProba) {
      bestRoot = root;
      max = rootProba;
    }
  }
  return bestRoot;
}

void AbstractReconciliationModel::updateRoot(pllmod_treeinfo_t &treeinfo) 
{
  setRoot(computeMLRoot(treeinfo));
}

double AbstractReconciliationModel::getSumLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  ScaledValue total;
  vector<pll_unode_t *> roots;
  getRoots(*treeinfo, roots, _geneIds);
  for (auto root: roots) {
    total += getRootLikelihood(*treeinfo, root);
  }
  return total.getLogValue(); 
}


void AbstractReconciliationModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo)
{
  vector<ScaledValue> zeros(speciesNodesCount_); 
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, _geneIds);
  for (auto root: roots) {
    pll_unode_t virtualRoot;
    virtualRoot.next = root;
    virtualRoot.node_index = root->node_index + _maxGeneId + 1;
    computeRootLikelihood(treeinfo, &virtualRoot);
  }
}
  
void AbstractReconciliationModel::inferMLScenario(shared_ptr<pllmod_treeinfo_t> treeinfo, Scenario &scenario)
{
  // make sure the CLVs are filled
  updateCLVs(*treeinfo);
  computeLikelihoods(*treeinfo); 
  
  pll_unode_t *geneRoot = 0;
  pll_rnode_t *speciesRoot = 0;
  computeMLRoot(*treeinfo, geneRoot, speciesRoot);
 /*
  pll_unode_t *geneRoot = computeMLRoot(*treeinfo);
  pll_rnode_t *speciesRoot = speciesTree_->root;
  */
  
  scenario.setGeneRoot(geneRoot);
  scenario.setSpeciesTree(speciesTree_);
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = geneRoot->node_index + _maxGeneId + 1;
  backtrace(&virtualRoot, speciesRoot, scenario, true);
}
  
