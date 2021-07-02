#include "PLLUnrootedTree.hpp"

#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>
#include <trees/PLLRootedTree.hpp>  
#include <stack>
#include <functional>
#include <sstream>


void defaultUnodePrinter(pll_unode_t *node, 
    std::stringstream &ss)
{
  if (node->label) {
    ss << node->label;
  }
  ss << ":" << node->length;
}


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

PLLUnrootedTree::PLLUnrootedTree(PLLRootedTree &rootedTree):
  _tree(pll_rtree_unroot(rootedTree.getRawPtr()), utreeDestroy)
{
  pll_unode_t *root = 
    _tree->nodes[_tree->tip_count + _tree->inner_count - 1];
  pll_utree_reset_template_indices(root, _tree->tip_count);
}
  
std::unique_ptr<PLLUnrootedTree> PLLUnrootedTree::buildFromStrOrFile(const std::string &strOrFile)
{
  std::unique_ptr<PLLUnrootedTree> res;
  try {
    res = std::make_unique<PLLUnrootedTree>(strOrFile, true);
  } catch (...) {
    try {
      res = std::make_unique<PLLUnrootedTree>(strOrFile, false);
    } catch (...) {
    }
  }
  return res;
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
  
CArrayRange<pll_unode_t*> PLLUnrootedTree::getLeaves() const
{
  return CArrayRange<pll_unode_t*>(_tree->nodes, getLeavesNumber());
}

CArrayRange<pll_unode_t*> PLLUnrootedTree::getNodes() const
{
  return CArrayRange<pll_unode_t*>(_tree->nodes, getNodesNumber());
}

CArrayRange<pll_unode_t*> PLLUnrootedTree::getInnerNodes() const
{
  return CArrayRange<pll_unode_t*>(_tree->nodes + getLeavesNumber(), getInnerNodesNumber());
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
  
pll_unode_t *PLLUnrootedTree::getNode(unsigned int node_index) const
{
  return _tree->nodes[node_index];
}

pll_unode_t *PLLUnrootedTree::getAnyInnerNode() const
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


static bool isBranchIn(pll_unode_t *b, 
    const std::unordered_set<pll_unode_t *> &branches)
{
  return branches.find(b) != branches.end() 
    || branches.find(b->back) != branches.end();
}

std::unordered_set<pll_unode_t *> PLLUnrootedTree::getBranches()
{
  std::unordered_set<pll_unode_t *> branches;
  for (auto node: getNodes()) {
    if (!isBranchIn(node, branches)) {
      branches.insert(node);
    }
    if (node->next) {
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
      node = node->next;
      if (!isBranchIn(node, branches)) {
        branches.insert(node);
      }
    }
  }
  return branches;
}


static void fillPostOrder(pll_unode_t *node,
    std::vector<pll_unode_t*> &nodes,
    std::vector<char> &markedNodes)
{
  // we already traversed this node
  if (markedNodes[node->node_index]) {
    return;
  }
  // mark the node as traversed
  markedNodes[node->node_index] = true;
  // first process children
  if (node->next) {
    fillPostOrder(node->next->back, nodes, markedNodes);
    fillPostOrder(node->next->next->back, nodes, markedNodes);
  }
  nodes.push_back(node);
}


std::vector<pll_unode_t*> PLLUnrootedTree::getPostOrderNodes(bool innerOnly)
{
  std::vector<pll_unode_t*> nodes;
  std::vector<char> markedNodes(getDirectedNodesNumber(), false);
  if (innerOnly) {
    for (auto node: getLeaves()) {
      markedNodes[node->node_index] = true;
    }
  }
  // do the post order traversal from all possible virtual roots 
  for (auto node: getLeaves()) {
    fillPostOrder(node->back, nodes, markedNodes);
  }
  if (innerOnly) {
    assert(nodes.size() == getDirectedNodesNumber() - getLeavesNumber());
  } else {
    assert(nodes.size() == getDirectedNodesNumber());
  }
  return nodes;
}

static void computePairwiseDistancesRec(pll_unode_t *currentNode, 
    double currentDistance,
    VectorDouble &distancesVector)
{
  currentDistance += currentNode->length;
  if (!currentNode->next) {
    // leaf
    distancesVector[currentNode->node_index] = currentDistance;
    return;
  }
  computePairwiseDistancesRec(currentNode->next->back, 
      currentDistance, 
      distancesVector);
  computePairwiseDistancesRec(currentNode->next->next->back, 
      currentDistance, 
      distancesVector);
} 

void PLLUnrootedTree::computePairwiseDistances(MatrixDouble &distances,
    bool leavesOnly)
{
  auto M = leavesOnly ? getLeavesNumber() : getDirectedNodesNumber(); 
  auto N = getLeavesNumber();
  auto MIter = leavesOnly ? getLeavesNumber() : getNodesNumber();
  VectorDouble zeros(N, 0.0);
  distances = MatrixDouble(M, zeros);
  for (unsigned int i = 0; i < MIter; ++i) {
    auto node = getNode(i);
    auto &distancesVector = distances[node->node_index];
    computePairwiseDistancesRec(node->back, 0.0, distancesVector);
    if (node->next) {
      // compute distances to leaves in all three directions
      computePairwiseDistancesRec(node->next->back, 
          0.0, 
          distancesVector);
      computePairwiseDistancesRec(node->next->next->back, 
          0.0, 
          distancesVector);
      // also update the two other directed nodes
      distances[node->next->node_index] = distancesVector;
      distances[node->next->next->node_index] = distancesVector;
    }
  }
}


static void getCladeRec(pll_unode_t *node, 
    std::unordered_set<unsigned int> &clade)
{
  if (node->next) {
    getCladeRec(node->next->back, clade);
    getCladeRec(node->next->next->back, clade);
  } else {
    clade.insert(node->node_index);
  }
}

std::unordered_set<unsigned int> 
  PLLUnrootedTree::getClade(pll_unode_t *node)
{
  std::unordered_set<unsigned int> clade;
  getCladeRec(node, clade);
  return clade;
}


// look for *pv (or one of his nexts) under the oriented node u.
// if it is found, *pv is updated with the found next, and 
// the function returns true (and false otherwise)
static bool orientAux(pll_unode_t *u, 
    pll_unode_t **pv,
    std::stack<pll_unode_t *> &path)
{
  path.push(u);
  auto v = *pv;
  if (v == u) {
    return true;
  }
  if (v->next) {
    if (v->next == u) {
      *pv = v->next;
      return true;
    } else if (v->next->next == u) {
      *pv = v->next->next;
      return true;
    }
  }
  if (!u->next) { 
    // end of recursion, we did not find *v
    path.pop();
    return false;
  } else {
    if (orientAux(u->next->back, pv, path) || 
      orientAux(u->next->next->back, pv, path)) {
      return true;
    } else {
      path.pop();
      return false;
    }
  }
}

static void stackToVector(std::stack<pll_unode_t *> s,
    std::vector<pll_unode_t *> &v)
{
  v.clear();
  v.resize(s.size());
  for (int i = s.size() - 1; i >= 0; --i) {
    v[i] = s.top();
    s.pop();
  }
}

void PLLUnrootedTree::orientTowardEachOther(pll_unode_t **pu,
    pll_unode_t **pv,
    std::vector<pll_unode_t *> &branchesPath)
{
  assert((*pu) != (*pv));
  auto *u = *pu;
  std::stack<pll_unode_t *> path;
  if (orientAux(u->back, pv, path)) {
    stackToVector(path, branchesPath);    
    return;
  }
  assert(path.size() == 0);
  if (u->next) {
    if (orientAux(u->next->back, pv, path)) {
      *pu = u->next;
      stackToVector(path, branchesPath);    
      return;
    } else if (orientAux(u->next->next->back, pv, path)) {
      *pu = u->next->next;
      stackToVector(path, branchesPath);    
      return;
    }
  }
  assert(false);
}



std::vector<double> PLLUnrootedTree::getMADRelativeDeviations()
{
  MatrixDouble distances;
  computePairwiseDistances(distances, false);
  auto nodes = getPostOrderNodes();  
  std::vector<double> deviations(nodes.size());
  for (auto node: getPostOrderNodes()) {
    // we reuse notations from the MAD paper (I, J, i, j, rho, b, c)
    // I and J are the set of leaves of each side of the branch node
    // rho is the relative position of the rooting that minimizes 
    // the squared relative deviation
    auto I = getClade(node); 
    auto J = getClade(node->back);
    auto rho = 0.0;
    auto rhoDen = 0.0;
    auto i = node->node_index;
    auto Dij = node->length;
    for (auto b: I) {
      for (auto c: J) {
        auto invPowBC = pow(distances[b][c], -2.0);
        rho += (distances[b][c] - 2.0 * distances[i][b]) * invPowBC;
        rhoDen += invPowBC;
      }
    }
    rho = rho / (2.0 * Dij * rhoDen);
    rho = std::max(std::min(rho, 1.0), 0.0);
    auto deviation = 0.0;
    for (auto b: I) {
      for (auto c: J) {
        auto v = (2.0 * (distances[i][b] + rho * Dij) / distances[b][c]);
        v -= 1.0;
        deviation += pow(v, 2.0);
      }
    }
    deviations[node->node_index] = deviation;
  }
  return deviations;
}

static void printAux(pll_unode_t *node,
    std::stringstream &ss,
    UnodePrinter f)
{
  if (node->next) {
    ss << "(";
    printAux(node->next->back, ss, f);
    ss << ",";
    printAux(node->next->next->back, ss, f);
    ss << ")";
  }
  f(node, ss);
}

std::string PLLUnrootedTree::getNewickString(UnodePrinter f,
      pll_unode_t *root, 
      bool rooted)
{
  std::stringstream ss;
  if (!root) {
    root = getAnyInnerNode();
  }
  if (rooted) {
    ss << "(";
    printAux(root, ss, f);
    ss << ",";
    printAux(root->back, ss, f);
    ss << ");";
  } else {
    ss << "(";
    printAux(root->back, ss, f);
    ss << ",";
    printAux(root->next->back, ss, f);
    ss << ",";
    printAux(root->next->next->back, ss, f);
    ss << ");";
  }
  return ss.str();
}

std::unordered_set<std::string> PLLUnrootedTree::getLabels() const
{
  std::unordered_set<std::string> res;
  for (auto node:  getLeaves()) {
    if (node->label) {
      res.insert(node->label);
    }
  }
  return res;
}

static pll_unode_t *searchForSet(pll_unode_t *node, 
    std::unordered_set<std::string> &currentNodeSet,
    const std::unordered_set<std::string> &set)
{
  if (!node->next) {
    if (node->label) {
      currentNodeSet.insert(std::string(node->label));
    }
  } else {
    auto left = node->next->back;
    auto right = node->next->next->back;
    std::unordered_set<std::string> rightNodeSet;
    auto res1 = searchForSet(left, currentNodeSet, set);
    if (res1) {
      return res1;
    }
    auto res2 = searchForSet(right, rightNodeSet, set);
    if (res2) {
      return res2;
    }
    for (auto &elem: rightNodeSet) {
      currentNodeSet.insert(elem);
    }
  }
  return (currentNodeSet == set) ? node : nullptr;
}


pll_unode_t *PLLUnrootedTree::getVirtualRoot(PLLRootedTree &referenceTree)
{
  std::unordered_set<std::string> leftRLeaves;
  PLLRootedTree::getLeafLabelsUnder(referenceTree.getRoot()->left, leftRLeaves);
  // find a leaf that is not in leftLeaves
  pll_unode_t *virtualRoot = nullptr;
  for (auto leaf: getLeaves()) {
    if (leftRLeaves.find(std::string(leaf->label)) == leftRLeaves.end()) {
      virtualRoot = leaf->back;
      break;
    }
  }
  std::unordered_set<std::string> temp;
  return searchForSet(virtualRoot, temp, leftRLeaves);
}

static size_t leafHash(pll_unode_t *leaf) {
  assert(leaf);
  std::hash<std::string> hash_fn;
  return hash_fn(std::string(leaf->label));
}

static size_t getTreeHashRec(pll_unode_t *node, size_t i) {
  assert(node);
  if (i == 0) 
    i = 1;
  if (!node->next) {
    return leafHash(node);
  }
  auto hash1 = getTreeHashRec(node->next->back, i + 1);
  auto hash2 = getTreeHashRec(node->next->next->back, i + 1);
  std::hash<size_t> hash_fn;
  auto m = std::min(hash1, hash2);
  auto M = std::max(hash1, hash2);
  return hash_fn(m * i + M);

}

static pll_unode_t *findMinimumHashLeafRec(pll_unode_t * root, size_t &hashValue)
{
  assert(root);
  if (!root->next) {
    hashValue = leafHash(root);
    return root;
  }
  auto n1 = root->next->back;
  auto n2 = root->next->next->back;
  size_t hash1, hash2;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    hashValue = hash1;
    return min1;
  } else {
    hashValue = hash2;
    return min2;
  }
}

static pll_unode_t *findMinimumHashLeaf(pll_unode_t * root) 
{
  assert(root);
  auto n1 = root;
  auto n2 = root->back;
  size_t hash1 = 0;
  size_t hash2 = 0;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    return min1;
  } else {
    return min2;
  }
}

size_t PLLUnrootedTree::getUnrootedTreeHash() const
{
  auto minHashLeaf = findMinimumHashLeaf(getAnyInnerNode());
  auto res = getTreeHashRec(minHashLeaf, 0) + getTreeHashRec(minHashLeaf->back, 0);
  return res;
}



static const char *orderChildren(pll_unode_t *node,
    std::vector<bool> &leftFirst)
{
  if (!node->next) {
    return node->label;
  }
  const char *label1 = orderChildren(node->next->back, leftFirst);
  const char *label2 = orderChildren(node->next->next->back, leftFirst);
  if (strcmp(label1, label2) < 0) {
    leftFirst[node->node_index] = true;
    return label1;
  } else {
    leftFirst[node->node_index] = false;
    return label2;
  }
}

static bool areIsomorphicAux(pll_unode_t *node1,
    pll_unode_t *node2,
    const std::vector<bool> &leftFirst1,
    const std::vector<bool> &leftFirst2)
{
  if (!node1->next || !node2->next) {
    // at least one is a leaf
    if (node1->next || node2->next) {
      // only one is a leaf 
      return false;
    }
    return strcmp(node1->label, node2->label) == 0;
  }
  // both are internal nodes
  auto l1 = node1->next->back;
  auto r1 = node1->next->next->back;
  auto l2 = node2->next->back;
  auto r2 = node2->next->next->back;
  if (!leftFirst1[node1->node_index]) {
    std::swap(l1, r1);
  }
  if (!leftFirst2[node2->node_index]) {
    std::swap(l2, r2);
  }
  return areIsomorphicAux(l1, l2, leftFirst1, leftFirst2) 
    && areIsomorphicAux(r1, r2, leftFirst1, leftFirst2);
}

bool PLLUnrootedTree::areIsomorphic(const PLLUnrootedTree &t1,
    const PLLUnrootedTree &t2)
{
  if (t1.getNodesNumber() != t2.getNodesNumber()) {
    return false;
  }
  auto startingLeaf1 = findMinimumHashLeaf(t1.getAnyInnerNode())->back;
  auto startingLeaf2 = findMinimumHashLeaf(t2.getAnyInnerNode())->back;
  std::vector<bool> leftFirst1(t1.getDirectedNodesNumber(), true);;
  std::vector<bool> leftFirst2(t2.getDirectedNodesNumber(), true);;
  orderChildren(startingLeaf1, leftFirst1);
  orderChildren(startingLeaf2, leftFirst2);
  return areIsomorphicAux(startingLeaf1, 
      startingLeaf2, 
      leftFirst1, 
      leftFirst2);
}

bool PLLUnrootedTree::isBinary() const
{
  for (auto node: getInnerNodes()) {
    assert(node->next);
    if (node->next->next->next != node) {
      return false;
    }
  }
  return true;
}
  
pll_unode_t *PLLUnrootedTree::findLeaf(const std::string &label)
{
  for (auto leaf: getLeaves()) {
    if (label == leaf->label) {
      return leaf;
    }
  }
  return nullptr;
}

