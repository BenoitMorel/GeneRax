#include <optimizers/DTLOptimizer.hpp>

#include <likelihoods/SubtreeRepeatsCache.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <limits>
#include <algorithm>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <iostream>
#include <cmath>

static bool isValidLikelihood(double ll) {
  return std::isnormal(ll) && ll < -0.0000001;
}


void updateLL(Parameters &rates, PerCoreGeneTrees &trees, pll_rtree_t *speciesTree, RecModel model) {
  rates.ensurePositivity();
  double ll = 0.0;
  for (auto &tree: trees.getTrees()) {
    ReconciliationEvaluation evaluation(speciesTree, tree.mapping, model, true);
    evaluation.setRates(rates);
    ll += evaluation.evaluate(*tree.geneTree);
  }
  ParallelContext::sumDouble(ll);
  if (!isValidLikelihood(ll)) {
    ll = -std::numeric_limits<double>::infinity();
  }
  rates.setScore(ll);
}

static bool lineSearchParameters(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, Parameters &currentRates, const Parameters &gradient, unsigned int &llComputationsLine)
{
  double alpha = 0.1; 
  double epsilon = 0.0000001;
  const double minAlpha = epsilon; //(dimensions == 2 ? epsilon : 0.001);
  const double minImprovement = 0.1;
  Parameters currentGradient(gradient);
  bool stop = true;
  while (alpha > minAlpha) {
    currentGradient.normalize(alpha);
    Parameters proposal = currentRates + (currentGradient * alpha);
    updateLL(proposal, geneTrees, speciesTree, model);
    llComputationsLine++;
    if (currentRates.getScore() + minImprovement < proposal.getScore()) {
      currentRates = proposal;
      stop = false;
      alpha *= 1.5;
    } else {
      alpha *= 0.5;
      if (!stop) {
        return !stop;
      }
    }
  }
  return !stop;
}

Parameters DTLOptimizer::optimizeParameters(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, const Parameters &startingParameters)
{
  double epsilon = 0.0000001;
  Parameters currentRates = startingParameters;
  updateLL(currentRates, geneTrees, speciesTree, model);
  unsigned int llComputationsGrad = 0;
  unsigned int llComputationsLine = 0;
  unsigned int dimensions = startingParameters.dimensions();
  Parameters gradient(dimensions);
  do {
    std::vector<Parameters> closeRates(dimensions, currentRates);
    for (unsigned int i = 0; i < dimensions; ++i) {
      Parameters closeRates = currentRates;
      closeRates[i] += epsilon;
      updateLL(closeRates, geneTrees, speciesTree, model);
      llComputationsGrad++;
      gradient[i] = (currentRates.getScore() - closeRates.getScore()) / (-epsilon);
    }
  } while (lineSearchParameters(geneTrees, speciesTree, model, currentRates, gradient, llComputationsLine));
  return currentRates;
}

Parameters DTLOptimizer::optimizeParametersGlobalDTL(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model)
{
  std::vector<Parameters> startingRates;
  if (Enums::freeParameters(model) == 2) {
    startingRates.push_back(Parameters(0.1, 0.2));
    startingRates.push_back(Parameters(0.2, 0.2));
    startingRates.push_back(Parameters(0.5, 0.5));
    startingRates.push_back(Parameters(0.5, 1.0));
    startingRates.push_back(Parameters(0.01, 0.01));
  } else {
    startingRates.push_back(Parameters(0.5, 0.5, 0.2));
    startingRates.push_back(Parameters(0.1, 0.2, 0.1));
    startingRates.push_back(Parameters(0.2, 0.2, 0.0));
    startingRates.push_back(Parameters(0.01, 0.01, 0.01));
  }
  Parameters best;
  best.setScore(-10000000000);
  for (auto rates: startingRates) {
    Parameters newRates = optimizeParameters(geneTrees, speciesTree, model, rates);
    bool stop = (fabs(newRates.getScore() - best.getScore()) < 3.0);
    if (newRates.getScore() > best.getScore()) {
      best = newRates;
    }
    if (stop) {
      break;
    }
  }
  return best; 
}


Parameters DTLOptimizer::optimizeParametersPerSpecies(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model)
{
  unsigned int speciesCount = speciesTree->tip_count + speciesTree->inner_count;
  Parameters globalRates = optimizeParametersGlobalDTL(geneTrees, speciesTree, model);
  Parameters startingSpeciesRates(speciesCount, globalRates);
  Parameters rates = DTLOptimizer::optimizeParameters(geneTrees, speciesTree, model, startingSpeciesRates);
  return rates; 
}




