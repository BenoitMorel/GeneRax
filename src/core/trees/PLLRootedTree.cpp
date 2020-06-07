#include "PLLRootedTree.hpp"

#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>
#include <set>
#include <cstring>
#include <maths/Random.hpp>


void PLLRootedTree::setSon(pll_rnode_t *parent, pll_rnode_t *newSon, bool left)
{
  newSon->parent = parent;
  if (left) {
    parent->left = newSon;
  } else {
    parent->right = newSon;
  }
}


static pll_rnode_t *createNode(const std::string &label, std::vector<pll_rnode_t *> &allNodes) {
  pll_rnode_t *node = static_cast<pll_rnode_t *>(malloc(sizeof(pll_rnode_t)));
  node->label = 0;
  if (label.size()) {
    node->label = static_cast<char *>(calloc(label.size() + 1, sizeof(char)));
    assert(node->label);
    strcpy(node->label, label.c_str());
  }
  node->node_index = static_cast<unsigned int>(allNodes.size());
  node->length = 0.1;
  node->parent = 0;
  node->left = 0;
  node->right = 0;
  node->data = 0;
  allNodes.push_back(node);
  return node;
}

static void destroyNodeData(void *)
{

}

void rtreeDestroy(pll_rtree_t *rtree) {
  if(!rtree)
    return;
  pll_rtree_destroy(rtree, destroyNodeData);
}

static pll_rtree_t *buildUtree(const std::string &str, bool isFile)
{
  if (isFile) {
    return LibpllParsers::readRootedFromFile(str);
  } else {
    return LibpllParsers::readRootedFromStr(str);
  }
}

PLLRootedTree::PLLRootedTree(const std::string &str, bool isFile):
  _tree(buildUtree(str, isFile), rtreeDestroy)
{
  setMissingLabels();
  setMissingBranchLengths();
}

PLLRootedTree::PLLRootedTree(const std::unordered_set<std::string> &labels):
  _tree(buildRandomTree(labels), rtreeDestroy)
{
  setMissingLabels();
  setMissingBranchLengths();
}

void PLLRootedTree::save(const std::string &fileName) const
{
  LibpllParsers::saveRtree(_tree->root, fileName);
}
  
std::string PLLRootedTree::getNewickString() const
{
  std::string res;
  LibpllParsers::getRtreeNewickString(_tree.get(), res);

  return res;
}
  
void PLLRootedTree::setMissingBranchLengths(double minBL)
{
  for (auto node: getNodes()) {
    if (0.0 == node->length) {
      node->length = minBL;
    } 
  }
}

void PLLRootedTree::setMissingLabels() 
{
  auto labels = getLabels(true);
  unsigned int i = 0;
  std::string prefix("s");
  for (auto node: getNodes()) {
    if (node->left) {
      std::string newLabel;
      if (node->label) {
        newLabel = std::string(node->label);
      }
      while (labels.find(newLabel) != labels.end() || newLabel.size() == 0) {
        newLabel = prefix + std::to_string(i++);
      }
      free(node->label);
      node->label = static_cast<char*>(malloc(sizeof(char) * (newLabel.size() + 1)));
      std::strcpy(node->label, newLabel.c_str());
    }
  }
}
  
CArrayRange<pll_rnode_t*> PLLRootedTree::getLeaves() const
{
  return CArrayRange<pll_rnode_t*>(_tree->nodes, getLeavesNumber());
}

CArrayRange<pll_rnode_t*> PLLRootedTree::getInnerNodes() const
{
  return CArrayRange<pll_rnode_t*>(_tree->nodes + getLeavesNumber(), getInnerNodesNumber());
}

CArrayRange<pll_rnode_t*> PLLRootedTree::getNodes() const
{
  return CArrayRange<pll_rnode_t*>(_tree->nodes, getNodesNumber());
}

unsigned int PLLRootedTree::getNodesNumber() const
{
  return getLeavesNumber() + getInnerNodesNumber();
}

unsigned int PLLRootedTree::getLeavesNumber() const
{
  return _tree->tip_count;
}

unsigned int PLLRootedTree::getInnerNodesNumber() const
{
  return _tree->inner_count;
}
  
pll_rnode_t *PLLRootedTree::getRoot() const
{
  return _tree->root;
}

pll_rnode_t *PLLRootedTree::getNode(unsigned int node_index) const
{
  return _tree->nodes[node_index];
}

pll_rnode_t *PLLRootedTree::getAnyInnerNode() const
{
  return getNode(getLeavesNumber());
}
  
std::unordered_set<std::string> PLLRootedTree::getLabels(bool leavesOnly) const
{
  std::unordered_set<std::string> res;
  for (auto node: (leavesOnly ? getLeaves() : getNodes())) {
    if (node->label) {
      res.insert(node->label);
    }
  }
  return res;
}

pll_rtree_t *PLLRootedTree::buildRandomTree(const std::unordered_set<std::string> &leafLabels)
{
  std::set<std::string> leaves;
  for (auto &leaf: leafLabels) {
    leaves.insert(leaf);
  }
  std::vector<pll_rnode_t *> allNodes;
  pll_rnode_t *root = 0;
  for (auto &label: leaves) {
    if (allNodes.size() == 0) {
      root = createNode(label, allNodes);
      continue;
    }
    auto brother = allNodes[static_cast<size_t>(Random::getInt()) % allNodes.size()];
    auto parent = createNode("", allNodes);
    auto node = createNode(label, allNodes);
    auto grandpa = brother->parent;
    if (grandpa) {
      setSon(grandpa, parent, grandpa->left == brother); 
    } else {
      root = parent;
    }
    bool randBool = static_cast<bool>(Random::getInt() % 2);
    setSon(parent, brother, randBool);
    setSon(parent, node, !randBool);
  }
  auto res = static_cast<pll_rtree_t *>(malloc(sizeof(pll_rtree_t)));
  res->root = root;
  res->nodes = static_cast<pll_rnode_t**>(malloc(sizeof(pll_rnode_t*) * allNodes.size()));
  for (unsigned int i = 0; i < allNodes.size(); ++i) {
    res->nodes[i] = allNodes[i];
  }
  res->tip_count = static_cast<unsigned int>(allNodes.size()) / 2 + 1;
  res->inner_count = static_cast<unsigned int>(allNodes.size()) / 2;
  res->edge_count = static_cast<unsigned int>(allNodes.size()) - 1;
  
  LibpllParsers::labelRootedTree(res);
  return res;
}

StringToUintMap PLLRootedTree::getLabelToIntMap()
{
  StringToUintMap map;
  for (auto node: getLeaves()) {
    map.insert({std::string(node->label), node->node_index});
  }
  return map;
}
 
static void fillPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t*> &nodes)
{
  if (node->left) {
    fillPostOrder(node->left, nodes);
    fillPostOrder(node->right, nodes);
  }
  nodes.push_back(node);
}

std::vector<pll_rnode_t*> PLLRootedTree::getPostOrderNodes() const
{ 
  std::vector<pll_rnode_t*> nodes;
  fillPostOrder(getRoot(), nodes);
  return nodes;
}



