#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <util/enums.hpp>
#include <maths/DTLRates.hpp>




class DTLOptimizer {
public:
  /**
   * Per-species rates optimization
   */
  static DTLRatesVector optimizeDTLRatesVector(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, DTLRatesVector *previousVector = 0);
  
  /**
   * Global rates optimization
   */
  static DTLRates optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model);
};


