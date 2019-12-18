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

  static Parameters optimizeParameters(PerCoreEvaluations &evaluations, 
      const Parameters &startingParameters);
  
  static Parameters optimizeParametersGlobalDTL(PerCoreEvaluations &evaluations, 
      Parameters *startingParameters = nullptr);

  static Parameters optimizeParametersPerSpecies(PerCoreEvaluations &evaluations, unsigned int speciesNodesNumber);

private:
  static void buildEvaluations(PerCoreGeneTrees &geneTrees, PLLRootedTree &speciesTree, RecModel model, Evaluations &evaluations);
};


