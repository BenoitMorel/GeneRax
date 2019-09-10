#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <util/enums.hpp>
#include <maths/Parameters.hpp>

class DTLOptimizer {
public:
  static Parameters optimizeParameters(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, const Parameters &startingParameters);
  static Parameters optimizeParametersGlobalDTL(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model);
  static Parameters optimizeParametersPerSpecies(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model);
};


