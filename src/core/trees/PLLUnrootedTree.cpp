#include "PLLUnrootedTree.hpp"

#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>

static void destroyNodeData(void *)
{

}

void utreeDestroy(pll_utree_t *utree) {
  if(!utree)
    return;
  pll_utree_destroy(utree, destroyNodeData);
}

static pll_utree_t *buildUtree(const std::string &str, bool isFile)
{
  if (isFile) {
    return LibpllParsers::readNewickFromFile(str);
  } else {
    return LibpllParsers::readNewickFromStr(str);
  }
}

PLLUnrootedTree::PLLUnrootedTree(const std::string &str, bool isFile):
  _tree(buildUtree(str, isFile), utreeDestroy)
{
}

PLLUnrootedTree::PLLUnrootedTree(const std::vector<const char*> &labels,
    unsigned int seed):
  _tree(pllmod_utree_create_random(static_cast<unsigned int>(labels.size()), &labels[0], seed), utreeDestroy)
{

}

void PLLUnrootedTree::save(const std::string &fileName)
{
  LibpllParsers::saveUtree(_tree->nodes[0], fileName, false);
}

void PLLUnrootedTree::setMissingBranchLengths(double minBL)
{
  for (auto node: getLeaves()) {
    if (0.0 == node->length) {
      node->length = minBL;
    } 
  }
  for (unsigned int i = _tree->tip_count; i < _tree->tip_count + _tree->inner_count; ++i) {
    if (0.0 == _tree->nodes[i]->length)
      _tree->nodes[i]->length = minBL;
    if (0.0 == _tree->nodes[i]->next->length)
      _tree->nodes[i]->next->length = minBL;
    if (0.0 == _tree->nodes[i]->next->next->length)
      _tree->nodes[i]->next->next->length = minBL;
  }  
}
  
CArrayRange<pll_unode_t*> PLLUnrootedTree::getLeaves()
{
  return CArrayRange<pll_unode_t*>(_tree->nodes, getLeavesNumber());
}

CArrayRange<pll_unode_t*> PLLUnrootedTree::getNodes()
{
  return CArrayRange<pll_unode_t*>(_tree->nodes, getNodesNumber());
}


unsigned int PLLUnrootedTree::getNodesNumber() const
{
  return getLeavesNumber() + getInnerNodesNumber();
}

unsigned int PLLUnrootedTree::getDirectedNodesNumber() const
{
  return getLeavesNumber() + 3 * getInnerNodesNumber();
}

unsigned int PLLUnrootedTree::getLeavesNumber() const
{
  return _tree->tip_count;
}

unsigned int PLLUnrootedTree::getInnerNodesNumber() const
{
  return _tree->inner_count;
}
  
pll_unode_t *PLLUnrootedTree::getNode(unsigned int node_index)
{
  return _tree->nodes[node_index];
}

pll_unode_t *PLLUnrootedTree::getAnyInnerNode()
{
  return getNode(getLeavesNumber());
}
  
std::unordered_set<std::string> PLLUnrootedTree::getLeavesLabels()
{
  std::unordered_set<std::string> res;
  for (auto leaf: getLeaves()) {
    if (leaf->label) {
      res.insert(std::string(leaf->label));
    }
  }
  return res;
}

static void fillPostOrder(pll_unode_t *node,
    std::vector<pll_unode_t*> &nodes,
    std::unordered_set<unsigned int> &markedNodes)
{
  // we already traversed this node
  if (markedNodes.find(node->node_index) != markedNodes.end()) {
    return;
  }
  // mark the node as traversed
  markedNodes.insert(node->node_index);
  // first process children
  if (node->next) {
    fillPostOrder(node->next->back, nodes, markedNodes);
    fillPostOrder(node->next->next->back, nodes, markedNodes);
  }
  nodes.push_back(node);
}

std::vector<pll_unode_t*> PLLUnrootedTree::getPostOrderNodes()
{
  std::vector<pll_unode_t*> nodes;
  std::unordered_set<unsigned int> markedNodes;
  // do the post order traversal from all possible virtual roots 
  for (auto node: getNodes()) {
    fillPostOrder(node, nodes, markedNodes);
    if (node->next) {
      fillPostOrder(node->next, nodes, markedNodes);
      fillPostOrder(node->next->next, nodes, markedNodes);
    }
  }
  assert(nodes.size() == getDirectedNodesNumber());
  return nodes;
}




