#include "SpeciesTree.hpp"
#include <cassert>

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <parallelization/ParallelContext.hpp>


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
  return getTaxaNumberAux(getRoot()) / 2;
}


void SpeciesTree::saveToFile(const std::string &newick)
{
  LibpllParsers::saveRtree(_speciesTree->root, newick);  
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

void setSon(pll_rnode_t *parent, pll_rnode_t *newSon, bool left)
{
  newSon->parent = parent;
  if (left) {
    parent->left = newSon;
  } else {
    parent->right = newSon;
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
  


