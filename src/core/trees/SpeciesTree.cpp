#include "SpeciesTree.hpp"
#include <cassert>

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <parallelization/ParallelContext.hpp>


SpeciesTree::SpeciesTree(const std::string &newick):
  _speciesTree(LibpllParsers::readRootedFromFile(newick))
{
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

bool SpeciesTree::canChangeRoot(bool left1, bool left2)
{
  assert(_speciesTree);
  auto root = _speciesTree->root;
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

void SpeciesTree::changeRoot(bool left1, bool left2)
{
  Logger::info << "Change root " << sideString(left1) << " " << sideString(left2) << std::endl;
  assert(canChangeRoot(left1, left2));
  auto root = _speciesTree->root;
  auto rootLeft = root->left;
  auto rootRight = root->right;
  auto A = rootLeft->left;
  auto B = rootLeft->right;
  auto C = rootRight->left;
  auto D = rootRight->right;

  if (left1 && left2) {
    _speciesTree->root = rootLeft;
    rootLeft->parent = 0;
    setSon(rootLeft, root, false);
    setSon(root, B, true);
    setSon(root, rootRight, false);
  } else if (!left1 && !left2) {
    _speciesTree->root = rootRight;
    rootRight->parent = 0;
    setSon(rootRight, root, true);
    setSon(root, C, false);
    setSon(root, rootLeft, true);
  } else if (left1 && !left2) {
    _speciesTree->root = rootLeft;
    rootLeft->parent = 0;
    setSon(rootLeft, rootLeft->right, true);
    setSon(rootLeft, root, false);
    setSon(root, A, false);
    setSon(root, rootRight, true);
  } else { // !left1 && left2
    _speciesTree->root = rootRight;
    rootRight->parent = 0;
    setSon(rootRight, root, true);
    setSon(rootRight, C, false);
    setSon(root, D, true);
    setSon(root, rootLeft, false);
  }


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
  return getTaxaNumberAux(getRoot());
}



