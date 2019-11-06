#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <util/enums.hpp>
#include <maths/Parameters.hpp>

class RootedTree;

class DTLOptimizer {
public:
  DTLOptimizer() = delete;

  static Parameters optimizeParameters(PerCoreGeneTrees &geneTrees, 
      PLLRootedTree &speciesTree, 
      RecModel model, 
      const Parameters &startingParameters);

  static Parameters optimizeParametersGlobalDTL(PerCoreGeneTrees &geneTrees, 
      PLLRootedTree &speciesTree, 
      RecModel model);

  static Parameters optimizeParametersPerSpecies(PerCoreGeneTrees &geneTrees, 
      PLLRootedTree &speciesTree, 
      RecModel model);
};


