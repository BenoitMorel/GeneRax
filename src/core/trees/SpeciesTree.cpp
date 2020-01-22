#include "SpeciesTree.hpp"
#include <cassert>

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/FileSystem.hpp>
#include <set>
#include <functional>


SpeciesTree::SpeciesTree(const std::string &newick, bool fromFile):
  _speciesTree(newick, fromFile)
{
}


SpeciesTree::SpeciesTree(const std::unordered_set<std::string> &leafLabels):
  _speciesTree(leafLabels)
{
}
  
std::unique_ptr<SpeciesTree> SpeciesTree::buildRandomTree() const
{
   return std::make_unique<SpeciesTree>(_speciesTree.getLabels(true));
}

  
SpeciesTree::SpeciesTree(const Families &families):
  _speciesTree(getLabelsFromFamilies(families))
{
}


std::string SpeciesTree::toString() const
{
  std::string newick;
  LibpllParsers::getRtreeHierarchicalString(_speciesTree.getRawPtr(), newick);
  return newick;
}


void SpeciesTree::saveToFile(const std::string &newick, bool masterRankOnly)
{
  if (masterRankOnly && ParallelContext::getRank()) {
    return;
  }
  _speciesTree.save(newick);
}
  
void SpeciesTree::addListener(Listener *listener)
{
  _listeners.push_back(listener);
}

void SpeciesTree::removeListener(Listener *listener)
{
  _listeners.erase(std::remove(_listeners.begin(), _listeners.end(), listener), _listeners.end());      
}

void SpeciesTree::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  for (auto listener: _listeners) {
    listener->onSpeciesTreeChange(nodesToInvalidate);
  }
}
  
static void setRootAux(SpeciesTree &speciesTree, pll_rnode_t *root) {
  speciesTree.getTree().getRawPtr()->root = root; 
  root->parent = 0;
}

bool SpeciesTreeOperator::canChangeRoot(const SpeciesTree &speciesTree, unsigned int direction)
{
  bool left1 = direction % 2;
  auto root = speciesTree.getTree().getRoot();
  assert(root);
  auto newRoot = left1 ? root->left : root->right;
  return newRoot->left && newRoot->right;
}

std::string sideString(bool left) {
  if (left) {
    return std::string("left");
  } else {
    return std::string("right");
  }
}


void SpeciesTreeOperator::changeRoot(SpeciesTree &speciesTree, unsigned int direction)
{
  bool left1 = direction % 2;
  bool left2 = direction / 2;
  assert(canChangeRoot(speciesTree, left1));
  auto root = speciesTree.getTree().getRoot();
  auto rootLeft = root->left;
  auto rootRight = root->right;
  auto A = rootLeft->left;
  auto B = rootLeft->right;
  auto C = rootRight->left;
  auto D = rootRight->right;
  std::unordered_set<pll_rnode_t *> nodesToInvalidate;
  nodesToInvalidate.insert(root);
  setRootAux(speciesTree, left1 ? rootLeft : rootRight);
  if (left1 && left2) {
    PLLRootedTree::setSon(rootLeft, root, false);
    PLLRootedTree::setSon(root, B, true);
    PLLRootedTree::setSon(root, rootRight, false);
  } else if (!left1 && !left2) {
    PLLRootedTree::setSon(rootRight, root, true);
    PLLRootedTree::setSon(root, C, false);
    PLLRootedTree::setSon(root, rootLeft, true);
  } else if (left1 && !left2) {
    PLLRootedTree::setSon(rootLeft, rootLeft->right, true);
    PLLRootedTree::setSon(rootLeft, root, false);
    PLLRootedTree::setSon(root, A, false);
    PLLRootedTree::setSon(root, rootRight, true);
  } else { // !left1 && left2
    PLLRootedTree::setSon(rootRight, root, true);
    PLLRootedTree::setSon(rootRight, C, false);
    PLLRootedTree::setSon(root, D, true);
    PLLRootedTree::setSon(root, rootLeft, false);
  }
  speciesTree.onSpeciesTreeChange(&nodesToInvalidate);
}

void SpeciesTreeOperator::revertChangeRoot(SpeciesTree &speciesTree, unsigned int direction)
{
  changeRoot(speciesTree, 3 - direction);
}

pll_rnode_t *getBrother(pll_rnode_t *node) {
  auto father = node->parent;
  assert(father);
  return father->left == node ? father->right : father->left;
}
  
bool SpeciesTreeOperator::canApplySPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft)
{
  auto pruneNode = speciesTree.getNode(prune);
  auto regraftNode = speciesTree.getNode(regraft);
  if (pruneNode->parent == regraftNode) {
    return false;
  }
  if (!pruneNode->parent || pruneNode == getBrother(pruneNode)) {
    return false;
  }
  while (regraftNode != speciesTree.getRoot()) {
    if (pruneNode == regraftNode) {
        return false;
    }
    regraftNode = regraftNode->parent;
  }
  return true;
}
  
unsigned int SpeciesTreeOperator::applySPRMove(SpeciesTree &speciesTree, 
    unsigned int prune, 
    unsigned int regraft)
{
  auto pruneNode = speciesTree.getNode(prune);
  auto pruneFatherNode = pruneNode->parent;
  assert(pruneFatherNode);
  auto pruneGrandFatherNode = pruneFatherNode->parent;
  auto pruneBrotherNode = getBrother(pruneNode);
  unsigned int res = pruneBrotherNode->node_index;
  std::unordered_set<pll_rnode_t *> nodesToInvalidate;
  // prune
  nodesToInvalidate.insert(pruneFatherNode);
  if (pruneGrandFatherNode) {
    nodesToInvalidate.insert(pruneGrandFatherNode);
    PLLRootedTree::setSon(pruneGrandFatherNode, pruneBrotherNode, pruneGrandFatherNode->left == pruneFatherNode);
  } else {
    setRootAux(speciesTree, pruneBrotherNode);
  }
  // regraft
  auto regraftNode = speciesTree.getNode(regraft);
  auto regraftParentNode = regraftNode->parent;
  if (!regraftParentNode) {
    // regraft is the root
    setRootAux(speciesTree, pruneFatherNode);
    PLLRootedTree::setSon(pruneFatherNode, regraftNode, pruneFatherNode->left != pruneNode);
  } else {
    PLLRootedTree::setSon(regraftParentNode, pruneFatherNode, regraftParentNode->left == regraftNode);
    PLLRootedTree::setSon(pruneFatherNode, regraftNode, pruneFatherNode->left != pruneNode);
    nodesToInvalidate.insert(regraftParentNode);
  }
  speciesTree.onSpeciesTreeChange(&nodesToInvalidate);
  return res;
}
  
void SpeciesTreeOperator::reverseSPRMove(SpeciesTree &speciesTree, 
    unsigned int prune, 
    unsigned int applySPRMoveReturnValue)
{
  applySPRMove(speciesTree, prune, applySPRMoveReturnValue);
}


// direction: 0 == from parent, 1 == from left, 2 == from right
static void recursiveGetNodes(pll_rnode_t *node, 
    unsigned int direction, 
    unsigned int radius, 
    std::vector<unsigned int> &nodes, 
    bool addNode = true)
{
  if (radius == 0 || node == 0) {
    return;
  }
  if (addNode) {
    nodes.push_back(node->node_index);
  }
  switch (direction) {
    case 0:
      recursiveGetNodes(node->left, 0, radius - 1, nodes);
      recursiveGetNodes(node->right, 0, radius - 1, nodes);
      break;
    case 1: case 2:
      recursiveGetNodes((direction == 1 ? node->right : node->left), 0, radius - 1, nodes);
      if (node->parent) {
        recursiveGetNodes(node->parent, node->parent->left == node ? 1 : 2, radius - 1, nodes);
      }
      break;
    default:
      assert(false);
  };

}
  
void SpeciesTreeOperator::getPossiblePrunes(SpeciesTree &speciesTree, std::vector<unsigned int> &prunes)
{
  for (auto pruneNode: speciesTree.getTree().getNodes()) {
    if (pruneNode == speciesTree.getTree().getRoot()) {
      continue;
    }
    prunes.push_back(pruneNode->node_index); 
  }
}
  
void SpeciesTreeOperator::getPossibleRegrafts(SpeciesTree &speciesTree, 
    unsigned int prune, 
    unsigned int radius, 
    std::vector<unsigned int> &regrafts)
{
  /**
   *  Hack: we do not add the nodes at the first radius, because they are equivalent to moves from the second radius
   */
  radius += 1;
  auto pruneNode = speciesTree.getNode(prune);
  auto pruneParentNode = pruneNode->parent;
  if (!pruneParentNode) {
    return;
  }
  if (pruneParentNode->parent) {
    auto parentDirection = static_cast<unsigned int>(pruneParentNode->parent->left == pruneParentNode ? 1 : 2);
    recursiveGetNodes(pruneParentNode->parent, parentDirection, radius, regrafts, false);
  }
  recursiveGetNodes(getBrother(pruneNode)->left, 0, radius, regrafts, false);
  recursiveGetNodes(getBrother(pruneNode)->right, 0, radius, regrafts, false);
}
  
static size_t leafHash(const pll_rnode_t *leaf) {
  assert(leaf);
  std::hash<std::string> hash_fn;
  return hash_fn(std::string(leaf->label));
}

static size_t getTreeHashRec(const pll_rnode_t *node, size_t i, bool useLeafHash) {
  assert(node);
  std::hash<size_t> hash_fn;
  if (i == 0) 
    i = 1;
  if (!node->left) {
    if (useLeafHash) {
      return leafHash(node);
    } else {
      return hash_fn(node->node_index);
    }
  }
  auto hash1 = getTreeHashRec(node->left, i + 1, useLeafHash);
  auto hash2 = getTreeHashRec(node->right, i + 1, useLeafHash);
  //Logger::info << "(" << hash1 << "," << hash2 << ") ";
  auto m = std::min(hash1, hash2);
  auto M = std::max(hash1, hash2);
  auto res = hash_fn(m * i + M);
  if (!useLeafHash) {
    res = hash_fn(res * i + node->node_index);
  }
  return res;
}
  
size_t SpeciesTree::getHash() const
{
  auto res = getTreeHashRec(getTree().getRoot(), 0, true);
  return res % 100000;  
}

size_t SpeciesTree::getNodeIndexHash() const
{
  auto res = getTreeHashRec(getTree().getRoot(), 0, false);
  return res % 100000;  
}


void SpeciesTree::getLabelsToId(std::unordered_map<std::string, unsigned int> &map) const
{
  map.clear();
  for (auto node: getTree().getNodes()) {
    map.insert(std::pair<std::string, unsigned int>(node->label, node->node_index));
  }
}

std::unordered_set<std::string> SpeciesTree::getLabelsFromFamilies(const Families &families)
{
  GeneSpeciesMapping mappings;
  for (const auto &family: families) {
    std::string geneTreeStr;
    FileSystem::getFileContent(family.startingGeneTree, geneTreeStr);
    mappings.fill(family.mappingFile, geneTreeStr);
  }
  std::unordered_set<std::string> leaves;
  for (auto &mapping: mappings.getMap()) {
    leaves.insert(mapping.second);
  }
  return leaves;
}
  

  
