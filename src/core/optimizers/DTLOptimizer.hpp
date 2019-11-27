#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <util/enums.hpp>
#include <maths/Parameters.hpp>
#include <memory>

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
      RecModel model,
      Parameters *startingParameters = nullptr);

  static Parameters optimizeParametersPerSpecies(PerCoreGeneTrees &geneTrees, 
      PLLRootedTree &speciesTree, 
      RecModel model);

private:
  static void buildEvaluations(PerCoreGeneTrees &geneTrees, PLLRootedTree &speciesTree, RecModel model, Evaluations &evaluations);
};


