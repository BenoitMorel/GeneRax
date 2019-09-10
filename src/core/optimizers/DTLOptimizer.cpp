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
    ll += evaluation.evaluate(tree.tree);
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

static Parameters optimizeParameters(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, const Parameters &startingParameters)
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
  
static void optimizeDTLRatesVectorGradient(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, DTLRatesVector &ratesVector)
{
  unsigned int dimensions = Enums::freeParameters(model);
  Parameters startingParameters(ratesVector.size() * dimensions);
  unsigned int index = 0;
  for (const auto &r: ratesVector.getRatesVector()) {
    for (unsigned int i = 0; i < dimensions; ++i) {
      startingParameters[index++] = r.rates[i];
    }
  }
  auto resultParameter = optimizeParameters(geneTrees, speciesTree, model, startingParameters);
  index = 0;
  for (auto &r: ratesVector.getRatesVector()) {
    for (unsigned int i = 0; i < dimensions; ++i) {
      r.rates[i] = resultParameter[index++];
    }
  }
  ratesVector.setLL(resultParameter.getScore());
}

DTLRatesVector DTLOptimizer::optimizeDTLRatesVector(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model)
{
  unsigned int speciesNumber = speciesTree->inner_count + speciesTree->tip_count;
  DTLRatesVector starting = DTLRatesVector(speciesNumber, optimizeDTLRates(geneTrees, speciesTree, model));
  optimizeDTLRatesVectorGradient(geneTrees, speciesTree, model, starting);
  return starting;
}



DTLRates DTLOptimizer::optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, const DTLRates &startingRates)
{

  unsigned int dimensions = Enums::freeParameters(model);
  DTLRates result;
  Parameters startingParameters(dimensions);
  for (unsigned int i = 0; i < dimensions; ++i) {
    startingParameters[i] = startingRates.rates[i];
  }
  auto resultParameter = optimizeParameters(geneTrees, speciesTree, model, startingParameters);
  for (unsigned int i = 0; i < dimensions; ++i) {
    result.rates[i] = resultParameter[i];
  }
  result.ll = resultParameter.getScore();
  return result;
}

double getRand() {
  return double(rand()) / double(RAND_MAX);
}

DTLRates DTLOptimizer::optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model)
{
  std::vector<DTLRates> startingRates;
  if (Enums::freeParameters(model) == 2) {
    startingRates.push_back(DTLRates(0.1, 0.2, 0.0));
    startingRates.push_back(DTLRates(0.2, 0.2, 0.0));
    startingRates.push_back(DTLRates(0.5, 0.5, 0.0));
    startingRates.push_back(DTLRates(0.5, 1.0, 0.0));
    startingRates.push_back(DTLRates(0.01, 0.01, 0.0));
  } else {
    startingRates.push_back(DTLRates(0.5, 0.5, 0.2));
    startingRates.push_back(DTLRates(0.1, 0.2, 0.1));
    startingRates.push_back(DTLRates(0.2, 0.2, 0.0));
    startingRates.push_back(DTLRates(0.01, 0.01, 0.01));
  }
  DTLRates best;
  best.ll = -10000000000;
  for (auto rates: startingRates) {
    DTLRates newRates = optimizeDTLRates(geneTrees, speciesTree, model, rates);
    bool stop = (fabs(newRates.ll - best.ll) < 3.0);
    if (newRates.ll > best.ll) {
      best = newRates;
    }
    if (stop) {
      break;
    }
  }
  return best; 
}
