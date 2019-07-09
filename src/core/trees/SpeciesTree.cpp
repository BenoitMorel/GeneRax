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
    ReconciliationEvaluation evaluation(_speciesTree, tree.mapping, model, true);
    evaluation.setRates(_rates.rates[0], _rates.rates[1], _rates.rates[2]);
    ll += evaluation.evaluate(tree.tree);
  }
  ParallelContext::sumDouble(ll);
  return ll;
}


