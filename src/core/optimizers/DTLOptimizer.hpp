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

struct OptimizationSettings {
  OptimizationSettings(): 
    lineSearchMinImprovement(0.1),
    optimizationMinImprovement(3.0),
    minAlpha(0.0000001),
    epsilon(0.0000001)
  {}



  double lineSearchMinImprovement;
  double optimizationMinImprovement;
  double minAlpha;
  double epsilon;

};

class DTLOptimizer {
public:
  DTLOptimizer() = delete;

  /**
   *  Generic parallel method for parameter optimization. Finds the 
   *  parameters (starting from startingParameters) that optimizes the 
   *  function equal to the sum (over parallel ranks) of the evaluations 
   *  functions. 
   *  @param evaluations the subset of functions allocated to the
   *                     current core
   *  @param startingParameters starting parameters
   *  @return The parameters that maximize the function
   */
  static Parameters optimizeParameters(PerCoreEvaluations &evaluations, 
      const Parameters &startingParameters,
      OptimizationSettings settings = OptimizationSettings());
 
  /**
   *  Finds the global parameters that maximize evaluations. Global 
   *  parameters means that they are not per-species nor per-families
   *  @param evaluations the subset of functions allocated to the 
   *                     current core
   *  @param startingParameters if not set, several preselected starting
   *                            parameters will be tried
   *  @return The parameters that maximize the function
   */
  static Parameters optimizeParametersGlobalDTL(PerCoreEvaluations &evaluations, 
      const Parameters *startingParameters = nullptr,
      OptimizationSettings settings = OptimizationSettings());

  /**
   * Same as optimizeParameters, but with a ModelParameters as input.
   */
  static ModelParameters optimizeModelParameters(PerCoreEvaluations &evaluations, 
      bool optimizeFromStartingParameters,
      const ModelParameters &startingParameters,
      OptimizationSettings settings = OptimizationSettings());

  /**
   * Finds the per-species parameters that maximize  evaluations
   *  @param evaluations the subset of functions allocated to the 
   *                     current core
   *  @param speciesNodesNumber number of species nodes
   *  @return The parameters that maximize the function
   */
  static Parameters optimizeParametersPerSpecies(PerCoreEvaluations &evaluations, unsigned int speciesNodesNumber);

private:
  static void buildEvaluations(PerCoreGeneTrees &geneTrees, PLLRootedTree &speciesTree, RecModel model, Evaluations &evaluations);
};


