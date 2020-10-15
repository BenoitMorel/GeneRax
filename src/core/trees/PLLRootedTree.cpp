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
  
pll_rnode_t *PLLRootedTree::getParent(unsigned int node_index) const 
{
  return getNode(node_index)->parent;
}

pll_rnode_t *PLLRootedTree::getNeighbor(unsigned int node_index) const
{
  auto node = getNode(node_index);
  auto parent = node->parent;
  assert(parent);
  return parent->left == node ? parent->right : parent->left;
}


pll_rnode_t *PLLRootedTree::getAnyInnerNode() const
{
  return getNode(getLeavesNumber());
}
  
void PLLRootedTree::getLeafLabelsUnder(pll_rnode_t *node,
    std::unordered_set<std::string> &labels)
{
  if (node->left) {
    getLeafLabelsUnder(node->left, labels);
    getLeafLabelsUnder(node->right, labels);
  } else {
    if (node->label) {
      labels.insert(std::string(node->label));
    }
  }
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
  
void PLLRootedTree::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *)
{
  if (_lcaCache) {
    buildLCACache();
  }
}
  
pll_rnode_t *PLLRootedTree::getLCA(pll_rnode_t *n1, pll_rnode_t *n2)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->lcas[n1->node_index][n2->node_index];
}
  
bool PLLRootedTree::areParents(pll_rnode_t *n1, pll_rnode_t *n2)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->parents[n1->node_index][n2->node_index];
}

std::vector<bool> &PLLRootedTree::getParentsCache(pll_rnode_t *n1)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->parents[n1->node_index];
}

std::vector<bool> &PLLRootedTree::getAncestorssCache(pll_rnode_t *n1)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->ancestors[n1->node_index];
}
  
static void fillWithChildren(pll_rnode_t *n1,
    pll_rnode_t *n2,
    std::vector<pll_rnode_t *> &n1lcas)
{
  if (!n2) {
    return;
  }
  n1lcas[n2->node_index] = n1;
  fillWithChildren(n1, n2->left, n1lcas);
  fillWithChildren(n1, n2->right, n1lcas);
}

/**
 * Recursion that starts with n2 == root and lca == root
 * traverse all nodes from the root to the leaves with n2
 * to fill n1 LCAs
 */
static void findn1LCAs(pll_rnode_t *n1,
    pll_rnode_t *n2,
    pll_rnode_t *lca,
    const std::unordered_set<pll_rnode_t*> &n1Ancestors,
    std::vector<pll_rnode_t *> &n1lcas)
{
  if (!n2) {
    // end of the recursion
    return;
  }
  if (n1 == n2) {
    // edge case: from now, the lca of n1 and all children of n2
    // is n1.
    fillWithChildren(n1, n2, n1lcas);
  } else {
    if (n1Ancestors.find(n2) != n1Ancestors.end()) {
      // n2 is an ancestor of n1, and thus the new lca
      // in the recursion
      lca = n2;
    }
    n1lcas[n2->node_index] = lca;
    findn1LCAs(n1, n2->left, lca, n1Ancestors, n1lcas);
    findn1LCAs(n1, n2->right, lca, n1Ancestors, n1lcas);
  }
}


static void findLCAs(pll_rnode_t *n1, std::vector<pll_rnode_t *> &n1lcas)
{
  std::unordered_set<pll_rnode_t *> n1Ancestors;
  auto it = n1;
  auto root = n1;
  while (it) {
    n1Ancestors.insert(it);
    root = it;
    it = it->parent;
  }
  findn1LCAs(n1, root, root, n1Ancestors, n1lcas);
}

void PLLRootedTree::buildLCACache()
{
  auto N = getNodesNumber();
  _lcaCache = std::make_unique<LCACache>();
  std::vector<pll_rnode_t *> nulls(N, nullptr);
  _lcaCache->lcas = std::vector<std::vector<pll_rnode_t *> >(N, nulls);
  std::vector<bool> falses(N, false);
  _lcaCache->parents = std::vector<std::vector<bool > >(N, falses);
  _lcaCache->ancestors = std::vector<std::vector<bool > >(N, falses);
  for (auto n: getNodes()) {
    findLCAs(n, _lcaCache->lcas[n->node_index]);
  }
  for (auto n1: getNodes()) {
    auto n2 = n1;
    while (n2) {
      _lcaCache->parents[n1->node_index][n2->node_index] = true;
      _lcaCache->parents[n2->node_index][n1->node_index] = true;
      _lcaCache->ancestors[n1->node_index][n2->node_index] = true;
      n2 = n2->parent;
    }
  }

}

StringToUint PLLRootedTree::getDeterministicLabelToId() const
{
  StringToUint res;
  auto v = getDeterministicIdToLabel();
  for (unsigned int i = 0; i < v.size(); ++i) {
    res[v[i]] = i;
  }
  return res;
}

std::vector<std::string> PLLRootedTree::getDeterministicIdToLabel() const
{
  std::vector<std::string> labels;
  for (auto node: getNodes()) {
    labels.push_back(std::string(node->label)); 
  }
  std::sort(labels.begin(), labels.end());
  return labels;
}
