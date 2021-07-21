#ifdef CODE_DISABLED
#include "SpeciesSplitScore.hpp"
#include <ccp/SpeciesSplitScore.hpp>

static const unsigned int INVALID_BID = static_cast<unsigned int>(-1);

static pll_rnode_t *getNeighbor(pll_rnode_t *node) {
  if (node->parent->left == node) {
    return node->parent->right;
  } else {
    return node->parent->left;
  }
}

/**
 *  Returns the children of node, assuming that 
 *  the root is at previousNode
 */ 
static void getChildren(pll_rnode_t *node,
    pll_rnode_t *previousNode,
    pll_rnode_t *& left,
    pll_rnode_t *& right)
{
  auto neighbor = getNeighbor(node);
  if (previousNode == node->parent || previousNode == neighbor) {
    // we go down the rooted tree
    left = node->left;
    right = node->right;
  } else if (node->parent->parent) {
    // we are not at the root
    // we go up and to the neighbor
    left = node->parent;
    right = neighbor;
  } else {
    // we reached the root, we go to the other side
    left = neighbor->left;
    right = neighbor->right;
  }
}
    
 
void SpeciesSplitScore::fillPathsRec(unsigned int fromSpid, 
    pll_rnode_t *node,
    pll_rnode_t* previousNode,
    BranchSet &path)
{
  if (previousNode && !node->left) { 
    // we reached the end of a path
    unsigned int toSpid = _labelToSpid[node->label];
    std::cout << "Path from " << fromSpid << " to " << node->label << ": " << path.count() << std::endl;
    _paths[fromSpid][toSpid] = path;
    return;
  }
  auto bid = _nodeIndexToBid[node->node_index];
  path.set(bid);
  pll_rnode_t *left;
  pll_rnode_t *right;
  getChildren(node, previousNode, left, right);
  fillPathsRec(fromSpid, left, node, path);
  fillPathsRec(fromSpid, right, node, path);
  path.unset(bid);
}


SpeciesSplitScore::SpeciesSplitScore(PLLRootedTree &speciesTree,
    const SpeciesSplits &splits):
  _nodeIndexToBid(speciesTree.getNodesNumber(), INVALID_BID),
  _emptyBranchSet(speciesTree.getLeavesNumber() - 3),
  _labelToSpid(splits.getLabelToSpid())
{
  auto leftRoot = speciesTree.getRoot()->left;
  auto rightRoot = speciesTree.getRoot()->right;
  auto topInternalBranch = leftRoot->left ? leftRoot : rightRoot;
  auto rootSonToDiscard = leftRoot->left ? rightRoot : leftRoot;
  for (auto node: speciesTree.getNodes()) {
    // map internal branches to branch ids
    // the root and it's right branch are excluded
    // because we represent an unrooted tree
    if (node == speciesTree.getRoot() 
        || node == speciesTree.getRoot()->left
        || node == rootSonToDiscard
        || !node->left) {
      continue;
    }
    unsigned int spid = _bidToNodeIndex.size();
    _bidToNodeIndex.push_back(node->node_index);
    _nodeIndexToBid[node->node_index] = spid;
  }
  if (rootSonToDiscard->left) {
    _nodeIndexToBid[rootSonToDiscard->node_index] = _nodeIndexToBid[topInternalBranch->node_index];
  }
  std::vector<BranchSet> emptySets(speciesTree.getLeavesNumber(), _emptyBranchSet);
  _paths.resize(speciesTree.getLeavesNumber());
  std::fill(_paths.begin(), _paths.end(), emptySets);
  _nodeIndexToBid[speciesTree.getRoot()->right->node_index] = 
    _nodeIndexToBid[speciesTree.getRoot()->left->node_index];
  assert(_bidToNodeIndex.size() == speciesTree.getLeavesNumber() - 3);
  for (auto leaf: speciesTree.getLeaves()) {
    BranchSet path = _emptyBranchSet;
    std::cout << "Starting from " << leaf->label << std::endl;
    fillPathsRec(_labelToSpid[leaf->label], leaf->parent, leaf, path);
  }
}



#endif
