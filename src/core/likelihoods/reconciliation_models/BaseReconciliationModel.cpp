#include "BaseReconciliationModel.hpp"

static bool fillNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd = nullptr)  
{
  bool addMyself = true;
  if (nodesToAdd) {
    addMyself = (nodesToAdd->find(node) != nodesToAdd->end());
  }
  if (node->left) {
    assert(node->right);
    addMyself |= fillNodesPostOrder(node->left, nodes, nodesToAdd);
    addMyself |= fillNodesPostOrder(node->right, nodes, nodesToAdd);
  }
  if (addMyself) {
    nodes.push_back(node);
  }
  return addMyself;
}


BaseReconciliationModel::BaseReconciliationModel(
    PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &recModelInfo):
  _info(recModelInfo),
  _speciesTree(speciesTree),
  _likelihoodMode(PartialLikelihoodMode::PartialGenes),
  _geneNameToSpeciesName(geneSpeciesMapping.getMap()),
  _allSpeciesNodesInvalid(true)
{
  initSpeciesTree();
  setFractionMissingGenes(_info.fractionMissingFile);
}


bool BaseReconciliationModel::fillPrunedNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd)
{
  bool addMyself = true;
  if (nodesToAdd) {
    addMyself = (nodesToAdd->find(node) != nodesToAdd->end());
  }
  if (getSpeciesLeft(node)) {
    assert(getSpeciesRight(node));
    addMyself |= fillPrunedNodesPostOrder(
        getSpeciesLeft(node), nodes, nodesToAdd);
    addMyself |= fillPrunedNodesPostOrder(
        getSpeciesRight(node), nodes, nodesToAdd);
  }
  if (addMyself) {
    nodes.push_back(node);
  }
  return addMyself;
}

void BaseReconciliationModel::initSpeciesTree()
{
  _allSpeciesNodesCount = _speciesTree.getNodesNumber();
  _speciesLeft = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesRight = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesParent = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesNameToId.clear();
  onSpeciesTreeChange(nullptr);
  for (auto node: _allSpeciesNodes) {
    if (!node->left) {
      _speciesNameToId[node->label] = node->node_index;
    }
  }
}

void BaseReconciliationModel::onSpeciesTreeChange(
    const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  if (!nodesToInvalidate) {
    _allSpeciesNodesInvalid = true;
  } else {
    assert(nodesToInvalidate->size());
    for (auto node: *nodesToInvalidate) {
      while (node) {
        _invalidatedSpeciesNodes.insert(node);
        node = node->parent;
      }
    }
    _invalidatedSpeciesNodes.insert(nodesToInvalidate->begin(), 
        nodesToInvalidate->end());
  }
  _allSpeciesNodes.clear();
  fillNodesPostOrder(_speciesTree.getRoot(), _allSpeciesNodes);
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    _speciesLeft[e] = speciesNode->left;
    _speciesRight[e] = speciesNode->right;
    _speciesParent[e] = speciesNode->parent;
  }
  _prunedRoot = _speciesTree.getRoot();
  if (_info.pruneSpeciesTree && _speciesCoverage.size()) {
    std::vector<pll_rnode_t *> pruned(_allSpeciesNodesCount, nullptr);
    for (auto speciesNode: _allSpeciesNodes) {
      auto e = speciesNode->node_index;
      if (!speciesNode->left) {
        pruned[e] = (_speciesCoverage[e] ? speciesNode : nullptr);   
      } else {
        auto left = _speciesLeft[e];
        auto right = _speciesRight[e];
        auto prunedLeft = pruned[left->node_index];
        auto prunedRight = pruned[right->node_index];
        if (prunedLeft && prunedRight) {
          _speciesLeft[e] = prunedLeft;
          _speciesRight[e] = prunedRight;
          pruned[e] = speciesNode;
          _speciesParent[prunedLeft->node_index] = speciesNode;
          _speciesParent[prunedRight->node_index] = speciesNode;
          _prunedRoot = speciesNode;
        } else if (!prunedLeft && prunedRight) {
          pruned[e] = prunedRight;
        } else if (prunedLeft && !prunedRight) {
          pruned[e] = prunedLeft;
        }
      }
    }
    _allSpeciesNodes.clear();
    fillPrunedNodesPostOrder(getPrunedRoot(), _allSpeciesNodes);
  }
  assert(_allSpeciesNodes.size()); // && _allSpeciesNodes.back() == _speciesTree.getRoot());
}

void BaseReconciliationModel::beforeComputeLogLikelihood()
{
  if (_allSpeciesNodesInvalid) { // update everything
    _speciesNodesToUpdate = _allSpeciesNodes;
  } else if (_invalidatedSpeciesNodes.size()) { // partial update
    // here, fill _speciesNodesToUpdate with the invalid nodes
    _speciesNodesToUpdate.clear();
    fillPrunedNodesPostOrder(getPrunedRoot(), _speciesNodesToUpdate, &_invalidatedSpeciesNodes);
  } else {
    _speciesNodesToUpdate.clear();
  }
  _allSpeciesNodesInvalid = false;
  _invalidatedSpeciesNodes.clear();
  //assert(!_speciesNodesToUpdate.size() || _speciesNodesToUpdate.back() == getPrunedRoot());
  recomputeSpeciesProbabilities();
}

void BaseReconciliationModel::setFractionMissingGenes(const std::string &fractionMissingFile)
{
  _fm = std::vector<double>(_allSpeciesNodesCount, 0.0);
  if (!fractionMissingFile.size()) {
    return;
  }
  std::ifstream is(fractionMissingFile);
  std::string species;
  double fm;
  while (is >> species >> fm) {
    _fm[_speciesNameToId[species]] = fm;
  }
}
  




