#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <util/enums.hpp>
#include <maths/DTLRates.hpp>
#include <maths/Parameters.hpp>



class DTLOptimizer {
public:
  /**
   * Per-species rates optimization
   */
  static DTLRatesVector optimizeDTLRatesVector(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model);
  
  /**
   * Global rates optimization
   */
  static DTLRates optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model);
  static DTLRates optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, const DTLRates &startingRates);
};


