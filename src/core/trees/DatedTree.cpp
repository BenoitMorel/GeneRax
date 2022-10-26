#include "DatedTree.hpp"
#include <iostream>
#include <IO/Logger.hpp>

DatedTree::DatedTree(PLLRootedTree *rootedTree, bool fromBL):
  _rootedTree(rootedTree)
{
  if (fromBL) {
    _orderedSpeciations = _rootedTree->getOrderedSpeciations();  
  } else {
    auto postOrder = rootedTree->getPostOrderNodes();
    for (auto it = postOrder.rbegin(); it != postOrder.rend(); ++it) {
      auto node = *it;
      if (node->left) {
        _orderedSpeciations.push_back(node);
      }
    }
  }
  _ranks.resize(rootedTree->getNodesNumber());
  unsigned int rank = 0;
  for (auto species: _orderedSpeciations) {
    _ranks[species->node_index] = rank++; 
  }
  for (auto leaf: this->_rootedTree->getLeaves()) {
    _ranks[leaf->node_index] = rank;
  }
}


void DatedTree::rescaleBranchLengths()
{
  assert(isConsistent());
  std::vector<double> heights(_rootedTree->getNodesNumber(), 0.0);
  double height = 0.0;
  
  for (auto node: _orderedSpeciations) {
    auto e = node->node_index;
    if (!node->parent) {
      node->length = 1.0;
      continue;
    }
    auto p = node->parent->node_index;
    node->length =  double(_ranks[e]) - double(_ranks[p]);
    height = _ranks[e];
  }
  height += 1.0;
  for (auto leaf: _rootedTree->getLeaves()) {
    auto p = leaf->parent->node_index;
    leaf->length = height - double(_ranks[p]);
  }

}

bool DatedTree::moveUp(unsigned int rank, bool force) 
{
  if (rank == 0) {
    return false;
  }
  return moveDown(rank - 1, force);
}

bool DatedTree::moveDown(unsigned int rank, bool force)
{
  assert(rank + 1 < _orderedSpeciations.size());
  auto n1 = _orderedSpeciations[rank];
  auto n2 = _orderedSpeciations[rank + 1];
  // n1 is above n2. We want to swap them
  if (!force && n2->parent == n1) {
    return false;
  }
  _orderedSpeciations[rank + 1] = n1;
  _orderedSpeciations[rank] = n2;
  assert(_ranks[n2->node_index] - _ranks[n1->node_index] == 1);
  _ranks[n1->node_index]++;
  _ranks[n2->node_index]--;
  return true;
}
  
void DatedTree::moveNodeToRoot(corax_rnode_t *node)
{
  auto e = node->node_index;
  auto rank = _ranks[e];
  while (rank != 0) {
    assert(moveUp(rank, true));
    
    rank = _ranks[e];
  }
  assert(rank == _ranks[e]);
}
 
bool DatedTree::isConsistent() const
{
  // check that _orderedSpeciations and ranks are consistent
  for (unsigned int i = 0; i < _orderedSpeciations.size() - 1; ++i) {
    auto n1 = _orderedSpeciations[i];
    auto n2 = _orderedSpeciations[i + 1];
    if (_ranks[n1->node_index] + 1 != _ranks[n2->node_index]) {
      return false;
    }
  }
  // check that the ranks are consistent with the tree structure
  for (auto node: _rootedTree->getNodes()) {
    if (node->parent) {
      auto p = node->parent->node_index;
      auto e = node->node_index;
      if (_ranks[p] > _ranks[e]) {
        return false;
      }
    }
  }
  return true;
}
  
void DatedTree::restore(const Backup &backup)
{
  _ranks = backup.ranks;
  auto speciations = _orderedSpeciations;
  for (auto node: speciations) {
    _orderedSpeciations[_ranks[node->node_index]] = node;
  }
}

