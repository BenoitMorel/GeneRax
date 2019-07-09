#pragma once

#include <IO/LibpllParsers.hpp>
#include <string>
#include <maths/DTLRates.hpp>
#include <util/enums.hpp>

class PerCoreGeneTrees;

class SpeciesTree {
public:
  SpeciesTree(const std::string &newick); 
  ~SpeciesTree();

  void setRates(const DTLRates &rates);
  const DTLRates &getRates() const;
  double computeReconciliationLikelihood(PerCoreGeneTrees &geneTrees, RecModel model);

private:
  pll_rtree_t *_speciesTree;
  DTLRates _rates;
};



