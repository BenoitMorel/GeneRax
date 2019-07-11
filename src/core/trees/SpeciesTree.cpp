#include "SpeciesTree.hpp"
#include <cassert>

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/FileSystem.hpp>

SpeciesTree::SpeciesTree(const std::string &newick, bool fromFile):
  _speciesTree(0)
{
  if (fromFile) {
    _speciesTree = LibpllParsers::readRootedFromFile(newick);
  } else {
    _speciesTree = LibpllParsers::readRootedFromStr(newick);
  }
  // make sure all nodes have labels
  LibpllParsers::labelRootedTree(_speciesTree);
  assert(_speciesTree);
}

void setSon(pll_rnode_t *parent, pll_rnode_t *newSon, bool left)
{
  newSon->parent = parent;
  if (left) {
    parent->left = newSon;
  } else {
    parent->right = newSon;
  }
}

pll_rnode_t *createNode(const std::string &label, std::vector<pll_rnode_t *> &allNodes) {
  pll_rnode_t *node = static_cast<pll_rnode_t *>(malloc(sizeof(pll_rnode_t)));
  node->label = 0;
  if (label.size()) {
    node->label = static_cast<char *>(calloc(label.size() + 1, sizeof(char)));
    assert(node->label);
    strcpy(node->label, label.c_str());
  }
  node->node_index = allNodes.size();
  node->length = 0.1;
  node->parent = 0;
  node->left = 0;
  node->right = 0;
  node->data = 0;
  allNodes.push_back(node);
  return node;
}

SpeciesTree::SpeciesTree(const std::unordered_set<std::string> &leafLabels)
{
  buildFromLabels(leafLabels);
}

void SpeciesTree::buildFromLabels(const std::unordered_set<std::string> &leafLabels)
{
  std::vector<pll_rnode_t *> allNodes;
  pll_rnode_t *root = 0;
  for (auto &label: leafLabels) {
    if (allNodes.size() == 0) {
      root = createNode(label, allNodes);
      continue;
    }
    auto brother = allNodes[rand() % allNodes.size()];
    auto parent = createNode("", allNodes);
    auto node = createNode(label, allNodes);
    auto grandpa = brother->parent;
    if (grandpa) {
      setSon(grandpa, parent, grandpa->left == brother); 
    } else {
      root = parent;
    }
    bool randBool = static_cast<bool>(rand() % 2);
    setSon(parent, brother, randBool);
    setSon(parent, node, !randBool);
  }
  _speciesTree = static_cast<pll_rtree_t *>(malloc(sizeof(pll_rtree_t)));
  _speciesTree->root = root;
  _speciesTree->nodes = static_cast<pll_rnode_t**>(malloc(sizeof(pll_rnode_t*) * allNodes.size()));
  for (unsigned int i = 0; i < allNodes.size(); ++i) {
    _speciesTree->nodes[i] = allNodes[i];
  }
  _speciesTree->tip_count = allNodes.size() / 2 + 1;
  _speciesTree->inner_count = allNodes.size() / 2;
  _speciesTree->edge_count = allNodes.size() - 1;
  
  assert(_speciesTree->tip_count == getTaxaNumber());
  LibpllParsers::labelRootedTree(_speciesTree);
}
  
SpeciesTree::SpeciesTree(const std::vector<FamiliesFileParser::FamilyInfo> &families)
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
  for (auto &label: leaves) {
    Logger::info << label << std::endl;
  }
  buildFromLabels(leaves);
}

SpeciesTree::~SpeciesTree()
{
  if (_speciesTree) {
    pll_rtree_destroy(_speciesTree, free);
  }
}

void SpeciesTree::setRates(const DTLRates &rates) 
{
  _rates = rates;
}

const DTLRates &SpeciesTree::getRates() const
{
  return _rates;
}
  

double SpeciesTree::computeReconciliationLikelihood(PerCoreGeneTrees &geneTrees, RecModel model)
{
  double ll = 0.0;
  for (auto &tree: geneTrees.getTrees()) {
    ReconciliationEvaluation evaluation(_speciesTree, tree.mapping, model, false);
    evaluation.setRates(_rates.rates[0], _rates.rates[1], _rates.rates[2]);
    ll += evaluation.evaluate(tree.tree);
  }
  ParallelContext::sumDouble(ll);
  return ll;
}

std::string SpeciesTree::toString() const
{
  std::string newick;
  LibpllParsers::getRtreeHierarchicalString(_speciesTree, newick);
  return newick;
}

unsigned int getTaxaNumberAux(const pll_rnode_t *node) 
{
  if (!node) {
    return 0;
  }
  return getTaxaNumberAux(node->left) + getTaxaNumberAux(node->right) + 1;
}

unsigned int SpeciesTree::getTaxaNumber() const
{
  return getTaxaNumberAux(getRoot()) / 2 + 1;
}


void SpeciesTree::saveToFile(const std::string &newick)
{
  LibpllParsers::saveRtree(_speciesTree->root, newick);  
}
  
pll_rnode_t *SpeciesTree::getRandomNode() 
{
  return getNode(rand() % (_speciesTree->tip_count + _speciesTree->inner_count)); 
}



bool SpeciesTreeOperator::canChangeRoot(const SpeciesTree &speciesTree, int direction)
{
  bool left1 = direction % 2;
  auto root = speciesTree.getRoot();
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


void SpeciesTreeOperator::changeRoot(SpeciesTree &speciesTree, int direction)
{
  bool left1 = direction % 2;
  bool left2 = direction / 2;
  assert(canChangeRoot(speciesTree, left1));
  auto root = speciesTree.getRoot();
  auto rootLeft = root->left;
  auto rootRight = root->right;
  auto A = rootLeft->left;
  auto B = rootLeft->right;
  auto C = rootRight->left;
  auto D = rootRight->right;
  speciesTree.setRoot(left1 ? rootLeft : rootRight);
  if (left1 && left2) {
    setSon(rootLeft, root, false);
    setSon(root, B, true);
    setSon(root, rootRight, false);
  } else if (!left1 && !left2) {
    setSon(rootRight, root, true);
    setSon(root, C, false);
    setSon(root, rootLeft, true);
  } else if (left1 && !left2) {
    setSon(rootLeft, rootLeft->right, true);
    setSon(rootLeft, root, false);
    setSon(root, A, false);
    setSon(root, rootRight, true);
  } else { // !left1 && left2
    setSon(rootRight, root, true);
    setSon(rootRight, C, false);
    setSon(root, D, true);
    setSon(root, rootLeft, false);
  }
}

void SpeciesTreeOperator::revertChangeRoot(SpeciesTree &speciesTree, int direction)
{
  changeRoot(speciesTree, 3 - direction);
}

pll_rnode_t *getBrother(pll_rnode_t *node) {
  auto father = node->parent;
  assert(father);
  return father->left == node ? father->right : father->left;
}
  
unsigned int SpeciesTreeOperator::applySPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft)
{
  auto pruneNode = speciesTree.getNode(prune);
  auto pruneFatherNode = pruneNode->parent;
  assert(pruneFatherNode);
  auto pruneGrandFatherNode = pruneFatherNode->parent;
  auto pruneBrotherNode = getBrother(pruneNode);
  unsigned int res = pruneBrotherNode->node_index;
  // prune
  if (pruneGrandFatherNode) {
    setSon(pruneGrandFatherNode, pruneBrotherNode, pruneGrandFatherNode->left == pruneFatherNode);
  } else {
    speciesTree.setRoot(pruneBrotherNode);
  }
  // regraft
  auto regraftNode = speciesTree.getNode(regraft);
  auto regraftParentNode = regraftNode->parent;
  if (!regraftParentNode) {
    // regraft is the root
    speciesTree.setRoot(pruneFatherNode);
    setSon(pruneFatherNode, regraftNode, pruneFatherNode->left != pruneNode);
  } else {
    setSon(regraftParentNode, pruneFatherNode, regraftParentNode->left == regraftNode);
    setSon(pruneFatherNode, regraftNode, pruneFatherNode->left != pruneNode);
  }
  return res;
}
  
void SpeciesTreeOperator::reverseSPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int applySPRMoveReturnValue)
{
  applySPRMove(speciesTree, prune, applySPRMoveReturnValue);
}


// direction: 0 == from parent, 1 == from left, 2 == from right
void recursiveGetNodes(pll_rnode_t *node, unsigned int direction, unsigned int radius, std::vector<unsigned int> &nodes)
{
  if (radius == 0 || node == 0) {
    return;
  }
  nodes.push_back(node->node_index);
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
  
void SpeciesTreeOperator::getPossibleRegrafts(SpeciesTree &speciesTree, unsigned int prune, unsigned int radius, std::vector<unsigned int> &regrafts)
{
  auto pruneNode = speciesTree.getNode(prune);
  auto pruneParentNode = pruneNode->parent;
  if (!pruneParentNode) {
    return;
  }
  if (pruneParentNode->parent) {
    int parentDirection = (pruneParentNode->parent->left == pruneParentNode ? 1 : 2);
    recursiveGetNodes(pruneParentNode->parent, parentDirection, radius, regrafts);
  }
  recursiveGetNodes(getBrother(pruneNode)->left, 0, radius, regrafts);
  recursiveGetNodes(getBrother(pruneNode)->right, 0, radius, regrafts);
}
  
void SpeciesTreeOptimizer::rootSlidingSearch(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model)
{
  double bestLL = speciesTree.computeReconciliationLikelihood(geneTrees, model);
  int bestMove = -1;
  do {
    bestMove = -1;
    Logger::info << "Current ll: " << bestLL << std::endl;
    for (unsigned int i = 0; i < 4; ++i) { 
      if (SpeciesTreeOperator::canChangeRoot(speciesTree, i)) {
        SpeciesTreeOperator::changeRoot(speciesTree, i); 
        double newLL = speciesTree.computeReconciliationLikelihood(geneTrees, model);
        Logger::info << "  New ll: " << newLL << std::endl;
        if (newLL > bestLL) {
          bestLL = newLL;
          bestMove = i;
        }
        SpeciesTreeOperator::revertChangeRoot(speciesTree, i); 
      }
    }
    if (bestMove != -1) {
      SpeciesTreeOperator::changeRoot(speciesTree, bestMove);  
    }
  } while (bestMove != -1); 
  Logger::info << "End of root sliding search: ll = " << bestLL << std::endl;
}
  
void rootExhaustiveSearchAux(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model, std::vector<unsigned int> &movesHistory, std::vector<unsigned int> &bestMovesHistory, double &bestLL, unsigned int &visits)
{
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  for (auto direction: moves) {
    if (SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      movesHistory.push_back(direction);
      SpeciesTreeOperator::changeRoot(speciesTree, direction);
      double ll = speciesTree.computeReconciliationLikelihood(geneTrees, model);
      visits++;
      if (ll > bestLL) {
        bestLL = ll;
        bestMovesHistory = movesHistory;
      }
      rootExhaustiveSearchAux(speciesTree, geneTrees, model, movesHistory, bestMovesHistory, bestLL, visits);
      SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
      movesHistory.pop_back();
    }
  }
}

void SpeciesTreeOptimizer::rootExhaustiveSearch(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model)
{
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  double bestLL = speciesTree.computeReconciliationLikelihood(geneTrees, model);
  unsigned int visits = 1;
  movesHistory.push_back(0);
  rootExhaustiveSearchAux(speciesTree, geneTrees, model, movesHistory, bestMovesHistory, bestLL, visits); 
  movesHistory[0] = 1;
  rootExhaustiveSearchAux(speciesTree, geneTrees, model, movesHistory, bestMovesHistory, bestLL, visits); 
  assert (visits == 2 * speciesTree.getTaxaNumber() - 3);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(speciesTree, bestMovesHistory[i]);
  }
}

