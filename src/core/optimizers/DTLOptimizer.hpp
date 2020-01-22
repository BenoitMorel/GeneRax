#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <util/enums.hpp>
#include <maths/Parameters.hpp>
#include <maths/ModelParameters.hpp>
#include <memory>

class RootedTree;

class DTLOptimizer {
public:
  DTLOptimizer() = delete;

  static Parameters optimizeParameters(PerCoreEvaluations &evaluations, 
      const Parameters &startingParameters);
  
  static Parameters optimizeParametersGlobalDTL(PerCoreEvaluations &evaluations, 
      const Parameters *startingParameters = nullptr);

  static ModelParameters optimizeModelParameters(PerCoreEvaluations &evaluations, 
      bool optimizeFromStartingParameters,
      const ModelParameters &startingParameters);

  static Parameters optimizeParametersPerSpecies(PerCoreEvaluations &evaluations, unsigned int speciesNodesNumber);

private:
  static void buildEvaluations(PerCoreGeneTrees &geneTrees, PLLRootedTree &speciesTree, RecModel model, Evaluations &evaluations);
};


