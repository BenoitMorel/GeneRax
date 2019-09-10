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
  

void updateLL(DTLRatesVector &rates, PerCoreGeneTrees &trees, std::vector<std::shared_ptr<ReconciliationEvaluation> > &evaluations) {
  rates.ensureValidity();
  rates.setLL(0.0);
  double ll = 0.0;
  for (unsigned int i = 0; i < trees.getTrees().size(); ++i) {
    auto &tree = trees.getTrees()[i];
    auto &evaluation = evaluations[i];
    evaluation->setRates(rates);
    ll += evaluation->evaluate(tree.tree);
  }
  ParallelContext::sumDouble(ll);
  if (!isValidLikelihood(ll)) {
    ll = -std::numeric_limits<double>::infinity();
  }
  rates.setLL(ll);
}


static bool lineSearchVector(PerCoreGeneTrees &geneTrees, std::vector<std::shared_ptr<ReconciliationEvaluation> > &evaluations, DTLRatesVector &currentRates, const DTLRatesVector &gradient, unsigned int &llComputationsLine)
{
  double alpha = 0.1; 
  double epsilon = 0.0000001;
  const double minAlpha = epsilon; //(dimensions == 2 ? epsilon : 0.001);
  const double minImprovement = 0.1;
  DTLRatesVector currentGradient(gradient);
  bool stop = true;
  while (alpha > minAlpha) {
    currentGradient.normalize(alpha);
    DTLRatesVector proposal = currentRates + (currentGradient * alpha);
    updateLL(proposal, geneTrees, evaluations);
    llComputationsLine++;
    if (currentRates.getLL() + minImprovement < proposal.getLL()) {
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
  Logger::info << resultParameter << std::endl;

  std::vector<std::shared_ptr<ReconciliationEvaluation> > evaluations;
  //SubtreeRepeatsCache cache;
  for (auto &tree: geneTrees.getTrees()) {
    auto evaluation = std::make_shared<ReconciliationEvaluation> (speciesTree, tree.mapping, model, true);
    evaluations.push_back(evaluation);
   // cache.addTree(tree.tree, tree.mapping);
   // evaluation->getReconciliationModel()->setCache(&cache); 
  }
  double epsilon = 0.000001;
  unsigned int species = ratesVector.size();
  DTLRatesVector currentRates = ratesVector;
  updateLL(currentRates, geneTrees, evaluations);
  unsigned int llComputationsLine = 0;
  DTLRatesVector gradient(species);
  do {
    for (unsigned int j = 0; j < species; ++j) {
      for (unsigned int i = 0; i < dimensions; ++i) {
        DTLRatesVector closeRates = currentRates;
        closeRates.getRates(j).rates[i] += epsilon;
        updateLL(closeRates, geneTrees, evaluations);
        gradient.getRates(j).rates[i] = (currentRates.getLL() - closeRates.getLL()) / (-epsilon);
      }
    }
  } while (lineSearchVector(geneTrees, evaluations, currentRates, gradient, llComputationsLine));
  //Logger::info << "Gradient ll: " << currentRates.getLL() << std::endl;
  ratesVector = currentRates;
  Logger::info << ratesVector << std::endl;
  
  for (auto &r: ratesVector.getRatesVector()) {
    for (unsigned int i = 0; i < dimensions; ++i) {
      r.rates[i] = startingParameters[index++];
    }
  }
  ratesVector.setLL(resultParameter.getScore());
  Logger::info << ratesVector << std::endl;
  ratesVector = currentRates;
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
  //Logger::info << "Best rates " << best << std::endl;
  return best; 
}
