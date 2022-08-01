#include "PLLRootedTree.hpp"

#include <IO/Logger.hpp>
#include <set>
#include <cstring>
#include <maths/Random.hpp>
#include <parallelization/ParallelContext.hpp>

static void * xmalloc(size_t size)
{ 
  void * t;
  t = malloc(size);
  return t;
} 
  
static char * xstrdup(const char * s)
{ 
  size_t len = strlen(s);
  char * p = (char *)xmalloc(len+1);
  return strcpy(p,s);
} 


static void dealloc_data(corax_rnode_t * node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

void corax_rtree_destroy(corax_rtree_t * tree,
                                  void (*cb_destroy)(void *))
{
  unsigned int i;
  corax_rnode_t * node;

  /* deallocate all nodes */
  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
  {
    node = tree->nodes[i];
    dealloc_data(node, cb_destroy);

    if (node->label)
      free(node->label);

    free(node);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}


static void fill_nodes_recursive(corax_unode_t * node,
                                 corax_unode_t ** array,
                                 unsigned int array_size,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index,
                                 unsigned int level)
{
  unsigned int index;
  if (!node->next)
  {
    /* tip node */
    index = *tip_index;
    *tip_index += 1;
  }
  else
  {
    /* inner node */
    corax_unode_t * snode = level ? node->next : node;
    do 
    {
      fill_nodes_recursive(snode->back, array, array_size, tip_index, 
                           inner_index, level+1);
      snode = snode->next;
    }
    while (snode != node);

    index = *inner_index;
    *inner_index += 1;
  }

  assert(index < array_size);
  array[index] = node;
}

static unsigned int utree_count_nodes_recursive(corax_unode_t * node, 
                                                unsigned int * tip_count,
                                                unsigned int * inner_count,
                                                unsigned int level)
{
  if (!node->next)
  {
    *tip_count += 1;
    return 1;
  }
  else
  {
    unsigned int count = 0;

    corax_unode_t * snode = level ? node->next : node;
	do 
	{
	  count += utree_count_nodes_recursive(snode->back, tip_count, inner_count, level+1);
	  snode = snode->next;
	}
	while (snode != node);

    *inner_count += 1;
	
	return count + 1;
  }
}

void PLLRootedTree::setSon(corax_rnode_t *parent, corax_rnode_t *newSon, bool left)
{
  newSon->parent = parent;
  if (left) {
    parent->left = newSon;
  } else {
    parent->right = newSon;
  }
}


static corax_rnode_t *createNode(const std::string &label, std::vector<corax_rnode_t *> &allNodes) {
  corax_rnode_t *node = static_cast<corax_rnode_t *>(malloc(sizeof(corax_rnode_t)));
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

void rtreeDestroy(corax_rtree_t *rtree) {
  if(!rtree)
    return;
  corax_rtree_destroy(rtree, destroyNodeData);
}

static corax_rtree_t *buildUtree(const std::string &str, bool isFile)
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
  ensureUniqueLabels();
  setMissingBranchLengths();
}

PLLRootedTree::PLLRootedTree(const std::unordered_set<std::string> &labels):
  _tree(buildRandomTree(labels), rtreeDestroy)
{
  ensureUniqueLabels();
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

void PLLRootedTree::ensureUniqueLabels() 
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
      labels.insert(newLabel);
    }
  }
}
  
CArrayRange<corax_rnode_t*> PLLRootedTree::getLeaves() const
{
  return CArrayRange<corax_rnode_t*>(_tree->nodes, getLeavesNumber());
}

CArrayRange<corax_rnode_t*> PLLRootedTree::getInnerNodes() const
{
  return CArrayRange<corax_rnode_t*>(_tree->nodes + getLeavesNumber(), getInnerNodesNumber());
}

CArrayRange<corax_rnode_t*> PLLRootedTree::getNodes() const
{
  return CArrayRange<corax_rnode_t*>(_tree->nodes, getNodesNumber());
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
  
corax_rnode_t *PLLRootedTree::getRoot() const
{
  return _tree->root;
}

corax_rnode_t *PLLRootedTree::getNode(unsigned int node_index) const
{
  return _tree->nodes[node_index];
}
  
corax_rnode_t *PLLRootedTree::getParent(unsigned int node_index) const 
{
  return getNode(node_index)->parent;
}

corax_rnode_t *PLLRootedTree::getNeighbor(unsigned int node_index) const
{
  auto node = getNode(node_index);
  auto parent = node->parent;
  assert(parent);
  return parent->left == node ? parent->right : parent->left;
}


corax_rnode_t *PLLRootedTree::getAnyInnerNode() const
{
  return getNode(getLeavesNumber());
}
  
void PLLRootedTree::getLeafLabelsUnder(corax_rnode_t *node,
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

corax_rtree_t *PLLRootedTree::buildRandomTree(const std::unordered_set<std::string> &leafLabels)
{
  std::set<std::string> leaves;
  for (auto &leaf: leafLabels) {
    leaves.insert(leaf);
  }
  std::vector<corax_rnode_t *> allNodes;
  corax_rnode_t *root = 0;
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
  auto res = static_cast<corax_rtree_t *>(malloc(sizeof(corax_rtree_t)));
  res->root = root;
  res->nodes = static_cast<corax_rnode_t**>(malloc(sizeof(corax_rnode_t*) * allNodes.size()));
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
 
static void fillPostOrder(corax_rnode_t *node, 
    std::vector<corax_rnode_t*> &nodes)
{
  if (node->left) {
    fillPostOrder(node->left, nodes);
    fillPostOrder(node->right, nodes);
  }
  nodes.push_back(node);
}

std::vector<corax_rnode_t*> PLLRootedTree::getPostOrderNodes() const
{ 
  std::vector<corax_rnode_t*> nodes;
  fillPostOrder(getRoot(), nodes);
  return nodes;
}
  
void PLLRootedTree::onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *)
{
  if (_lcaCache) {
    buildLCACache();
  }
}
  
corax_rnode_t *PLLRootedTree::getLCA(corax_rnode_t *n1, corax_rnode_t *n2)
{
  return getLCA(n1->node_index, n2->node_index);
}
  
corax_rnode_t *PLLRootedTree::getLCA(unsigned int nodeIndex1, unsigned int nodeIndex2)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->lcas[nodeIndex1][nodeIndex2];
  
}
  
bool PLLRootedTree::areParents(corax_rnode_t *n1, corax_rnode_t *n2)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->parents[n1->node_index][n2->node_index];
}

std::vector<bool> &PLLRootedTree::getParentsCache(corax_rnode_t *n1)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->parents[n1->node_index];
}

std::vector<bool> &PLLRootedTree::getAncestorssCache(corax_rnode_t *n1)
{
  if (!_lcaCache) {
    buildLCACache();
  }
  return _lcaCache->ancestors[n1->node_index];
}
  
static void fillWithChildren(corax_rnode_t *n1,
    corax_rnode_t *n2,
    std::vector<corax_rnode_t *> &n1lcas)
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
static void findn1LCAs(corax_rnode_t *n1,
    corax_rnode_t *n2,
    corax_rnode_t *lca,
    const std::unordered_set<corax_rnode_t*> &n1Ancestors,
    std::vector<corax_rnode_t *> &n1lcas)
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


static void findLCAs(corax_rnode_t *n1, std::vector<corax_rnode_t *> &n1lcas)
{
  std::unordered_set<corax_rnode_t *> n1Ancestors;
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
  std::vector<corax_rnode_t *> nulls(N, nullptr);
  _lcaCache->lcas = std::vector<std::vector<corax_rnode_t *> >(N, nulls);
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
  
std::vector<unsigned int> 
PLLRootedTree::getNodeIndexMapping(PLLRootedTree &otherTree)
{
  assert(otherTree.getNodesNumber() == getNodesNumber());
  std::map<std::string, corax_rnode_t*> otherTreeLabelToNode;
  for (auto node: otherTree.getLeaves()) {
    std::string label(node->label);
    otherTreeLabelToNode[label] = node;
  }
  std::vector<unsigned int> mapping(getNodesNumber(), 0);
  for (auto node: this->getPostOrderNodes()) {
    corax_rnode_t *otherNode = nullptr;
    if (!node->left) {
      // leaf case
      std::string label(node->label);
      otherNode = otherTreeLabelToNode[label];
    } else {
      auto leftId = node->left->node_index;
      auto rightId = node->right->node_index;
      auto otherLeft = otherTree.getNode(mapping[leftId]);
      auto otherRight = otherTree.getNode(mapping[rightId]);
      assert(otherLeft->parent == otherRight->parent);
      otherNode = otherLeft->parent;
    }
    mapping[node->node_index] = otherNode->node_index;
  }
  return mapping;
}

static corax_unode_t * rtree_unroot(corax_rnode_t * root, corax_unode_t * back)
{
  corax_unode_t * uroot = (corax_unode_t *)calloc(1,sizeof(corax_unode_t));
  uroot->back = back;
  uroot->label = (root->label) ? xstrdup(root->label) : NULL;
  uroot->length = uroot->back->length;

  if (!root->left)
  {
    uroot->next = NULL;
    return uroot;
  }

  uroot->next = (corax_unode_t *)calloc(1,sizeof(corax_unode_t));

  uroot->next->next = (corax_unode_t *)calloc(1,sizeof(corax_unode_t));
  uroot->next->next->next = uroot;

  uroot->next->length = root->left->length;
  uroot->next->back = rtree_unroot(root->left, uroot->next);
  uroot->next->next->length = root->right->length;
  uroot->next->next->back = rtree_unroot(root->right, uroot->next->next);

  return uroot;
}

corax_utree_t * corax_rtree_unroot(corax_rtree_t * tree)
{
  corax_rnode_t * root = tree->root;

  corax_rnode_t * new_root;

  corax_unode_t * uroot = (corax_unode_t*)calloc(1,sizeof(corax_unode_t));

  uroot->next = (corax_unode_t *)calloc(1,sizeof(corax_unode_t));

  uroot->next->next = (corax_unode_t *)calloc(1,sizeof(corax_unode_t));

  uroot->next->next->next = uroot;
  uroot->length = root->left->length + root->right->length;

  /* get the first root child that has descendants and make  it the new root */
  if (root->left->left)
  {
    new_root = root->left;
    uroot->back = rtree_unroot(root->right,uroot);
    /* TODO: Need to clean uroot in case of error */
    if (!uroot->back) return NULL;
  }
  else
  {
    new_root = root->right;
    uroot->back = rtree_unroot(root->left,uroot);
    /* TODO: Need to clean uroot in case of error*/
    if (!uroot->back) return NULL;
  }

  uroot->label = (new_root->label) ? xstrdup(new_root->label) : NULL;

  uroot->next->label = uroot->label;
  uroot->next->length = new_root->left->length;
  uroot->next->back = rtree_unroot(new_root->left, uroot->next);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->back) return NULL;

  uroot->next->next->label = uroot->label;
  uroot->next->next->length = new_root->right->length;
  uroot->next->next->back = rtree_unroot(new_root->right, uroot->next->next);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->next->back) return NULL;

  return corax_utree_wraptree(uroot,0);
}
  
  
 using LabelNodeIndex = std::pair<std::string, unsigned int>;

bool PLLRootedTree::areNodeIndicesParallelConsistent() const
{
  bool ok = true;
  std::vector<LabelNodeIndex> toSort;
  for (auto node: getPostOrderNodes()) {
    toSort.push_back(LabelNodeIndex(node->label, node->node_index));
    ok &= ParallelContext::isIntEqual(node->node_index);
  }
  std::sort(toSort.begin(), toSort.end());
  for (auto t: toSort) {
    ok &= ParallelContext::isIntEqual(t.second);
  }
  return ok;
}

