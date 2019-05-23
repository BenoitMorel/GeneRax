#include "AbstractReconciliationModel.hpp"
#include <IO/Logger.hpp>



  
AbstractReconciliationModel::AbstractReconciliationModel():
  rootedGeneTree_(false),
  geneRoot_(0),
  firstCall_(true),
  _maxGeneId(1)
{
}
void AbstractReconciliationModel::init(pll_rtree_t *speciesTree, const GeneSpeciesMapping &geneSpeciesMapping, bool rootedGeneTree)
{
  rootedGeneTree_ = rootedGeneTree;
  setSpeciesTree(speciesTree);
  geneNameToSpeciesName_ = geneSpeciesMapping.getMap();
}

void AbstractReconciliationModel::initFromUtree(pll_utree_t *tree) {
  auto treeSize = tree->tip_count + tree->inner_count;
  auto nodesNumber = tree->tip_count + 3 * tree->inner_count;
  _geneIds.clear();
  _allNodes.resize(nodesNumber);
  std::vector<bool> marked(nodesNumber, false);
  for (unsigned int i = 0; i < treeSize; ++i) {
    auto node = tree->nodes[i];
    _allNodes[node->node_index] = node;
    _geneIds.push_back(node->node_index);
    if  (node->next) {
      node = node->next;
      _allNodes[node->node_index] = node;
      _geneIds.push_back(node->node_index);
      node = node->next;
      _allNodes[node->node_index] = node;
      _geneIds.push_back(node->node_index);
    }
  }
}


void AbstractReconciliationModel::mapGenesToSpecies()
{
  geneToSpecies_.resize(_allNodes.size());
  for (auto node: _allNodes) {
    if (!node->next) {
      std::string speciesName = geneNameToSpeciesName_[std::string(node->label)]; 
      geneToSpecies_[node->node_index] = speciesNameToId_[speciesName];
    }
  }
  _cache.setGenesToSpecies(geneToSpecies_);
}

void AbstractReconciliationModel::setInitialGeneTree(pll_utree_t *tree)
{
  initFromUtree(tree);
  mapGenesToSpecies();
  _maxGeneId = static_cast<unsigned int>(_allNodes.size() - 1);
  invalidateAllCLVs();
  _cache.resetCache();
}

void AbstractReconciliationModel::fillNodesPostOrder(pll_rnode_t *node, std::vector<pll_rnode_t *> &nodes) 
{
  if (node->left) {
    assert(node->right);
    fillNodesPostOrder(node->left, nodes);
    fillNodesPostOrder(node->right, nodes);
  }
  nodes.push_back(node);
}

void AbstractReconciliationModel::setRates(double dupRate, 
  double lossRate,
  double transferRate)
{
  std::vector<double> dupRates(speciesNodesCount_, dupRate);
  std::vector<double> lossRates(speciesNodesCount_, lossRate);
  std::vector<double> transferRates(speciesNodesCount_, transferRate);
  setRates(dupRates, lossRates, transferRates);
  _cache.resetCache();
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
  _cache.resetCache(); 
}

void AbstractReconciliationModel::getRoots(std::vector<pll_unode_t *> &roots,
    const std::vector<unsigned int> &geneIds)
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
  std::vector<bool> marked(geneIds.size(), false);
  for (auto id: geneIds) {
    auto node = _allNodes[id];
    if (marked[node->node_index] || marked[node->back->node_index]) {
      continue;
    }
    roots.push_back(node->back);
    marked[node->node_index] = true;
  }
}
  
double AbstractReconciliationModel::computeLogLikelihood(pll_utree_t *tree)
{
  if (firstCall_) {
    setInitialGeneTree(tree);
    firstCall_ = false;
  }
  auto root = getRoot();
  updateCLVs();
  computeLikelihoods();
  if (rootedGeneTree_) {
    setRoot(computeMLRoot());
    while (root != getRoot()) {
      updateCLVs();
      computeLikelihoods();
      root = getRoot();
      setRoot(computeMLRoot());
    }
  }
  return getSumLikelihood();
}

pll_unode_t *AbstractReconciliationModel::getLeft(pll_unode_t *node, bool virtualRoot) const
{
  return virtualRoot ? node->next : node->next->back;
}

pll_unode_t *AbstractReconciliationModel::getRight(pll_unode_t *node, bool virtualRoot) const
{
  return virtualRoot ? node->next->back : node->next->next->back;
}

pll_unode_t *AbstractReconciliationModel::getLeftRepeats(pll_unode_t *node, bool virtualRoot)
{
  /*
  Logger::info << "getLeft " << node->node_index << " ";
  if (node->next) {
    Logger::info << node->next->back->node_index << " " << node->next->next->back->node_index;
  }
  Logger::info << std::endl;
  */
  return _cache.getRepeat(getLeft(node, virtualRoot));
}

pll_unode_t *AbstractReconciliationModel::getRightRepeats(pll_unode_t *node, bool virtualRoot)
{
  /*
  Logger::info << "getRight " << node->node_index << " ";
  if (node->next) {
    Logger::info << node->next->back->node_index << " " << node->next->next->back->node_index;
  }
  Logger::info << std::endl;
  */
  return _cache.getRepeat(getRight(node, virtualRoot));
}

void AbstractReconciliationModel::markInvalidatedNodesRec(pll_unode_t *node)
{
  _isCLVUpdated[node->node_index] = false;
  if (node->back->next) {
    markInvalidatedNodesRec(node->back->next);
    markInvalidatedNodesRec(node->back->next->next);
  }
}

void AbstractReconciliationModel::markInvalidatedNodes()
{
  for (auto nodeIndex: _invalidatedNodes) {
    auto node = _allNodes[nodeIndex];
    markInvalidatedNodesRec(node);
  }
  _invalidatedNodes.clear();
}

void AbstractReconciliationModel::updateCLVsRec(pll_unode_t *node)
{
  if (_isCLVUpdated[node->node_index]) {
    return;
  }
  std::stack<pll_unode_t *> nodes;
  nodes.push(node);
  while (!nodes.empty()) {
    auto currentNode = nodes.top();
    if (currentNode->next) {
      bool waitForChildren = false;
      auto left = getLeft(currentNode, false);
      auto right = getRight(currentNode, false);
      if (!_isCLVUpdated[left->node_index]) {
        nodes.push(left);
        waitForChildren = true;
      }
      if (!_isCLVUpdated[right->node_index]) {
        nodes.push(right);
        waitForChildren = true;
      }
      if (waitForChildren) {
        continue;
      }
    }
    updateCLV(currentNode);
    nodes.pop();
    _isCLVUpdated[currentNode->node_index] = true;
  }
}

void AbstractReconciliationModel::updateCLVs()
{
  markInvalidatedNodes();
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  for (auto root: roots) {
    updateCLVsRec(root);
    updateCLVsRec(root->back);
  }
}

void AbstractReconciliationModel::invalidateCLV(unsigned int nodeIndex)
{
  _invalidatedNodes.insert(nodeIndex);
}
  
void AbstractReconciliationModel::invalidateAllCLVs()
{
  _cache.resetCache();
  _isCLVUpdated = std::vector<bool>(_maxGeneId + 1, false);
}

void AbstractReconciliationModel::computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot) 
{
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  ScaledValue max;
  for (auto root: roots) {
    for (auto speciesNode: speciesNodes_) {
      ScaledValue ll = getRootLikelihood(root, speciesNode);
      if (max < ll) {
        max = ll;
        bestGeneRoot = root;
        bestSpeciesRoot = speciesNode;
      }
    }
  }
}

pll_unode_t *AbstractReconciliationModel::computeMLRoot()
{
  pll_unode_t *bestRoot = 0;
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  ScaledValue max;
  for (auto root: roots) {
    ScaledValue rootProba = getRootLikelihood(root);
    if (max < rootProba) {
      bestRoot = root;
      max = rootProba;
    }
  }
  return bestRoot;
}

double AbstractReconciliationModel::getSumLikelihood()
{
  ScaledValue total;
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  for (auto root: roots) {
    total += getRootLikelihood(root);
  }
  return total.getLogValue(); 
}


void AbstractReconciliationModel::computeLikelihoods()
{
  std::vector<ScaledValue> zeros(speciesNodesCount_); 
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  for (auto root: roots) {
    pll_unode_t virtualRoot;
    virtualRoot.next = root;
    virtualRoot.node_index = root->node_index + _maxGeneId + 1;
    computeRootLikelihood(&virtualRoot);
  }
}
  
void AbstractReconciliationModel::inferMLScenario(Scenario &scenario)
{
  // make sure the CLVs are filled
  updateCLVs();
  computeLikelihoods(); 
  
  pll_unode_t *geneRoot = 0;
  pll_rnode_t *speciesRoot = 0;
  computeMLRoot(geneRoot, speciesRoot);
  scenario.setGeneRoot(geneRoot);
  scenario.setSpeciesTree(speciesTree_);
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = geneRoot->node_index + _maxGeneId + 1;
  backtrace(&virtualRoot, speciesRoot, scenario, true);
}
  
