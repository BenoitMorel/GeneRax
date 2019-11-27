#include <optimizers/DTLOptimizer.hpp>

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


static void updateLL(Parameters &rates, Evaluations &evaluations) {
  rates.ensurePositivity();
  double ll = 0.0;
  for (auto evaluation: evaluations) {
    evaluation->setRates(rates);
    ll += evaluation->evaluate();
  }
  ParallelContext::sumDouble(ll);
  if (!isValidLikelihood(ll)) {
    ll = -std::numeric_limits<double>::infinity();
  }
  rates.setScore(ll);
}

static bool lineSearchParameters(Evaluations &evaluations, 
    Parameters &currentRates, 
    const Parameters &gradient, 
    unsigned int &llComputationsLine)
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
    updateLL(proposal, evaluations);
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

Parameters DTLOptimizer::optimizeParameters(PerCoreGeneTrees &geneTrees,
    PLLRootedTree &speciesTree, 
    RecModel model, const Parameters &startingParameters)
{
  Evaluations evaluations;
  buildEvaluations(geneTrees, speciesTree, model, evaluations); 
  double epsilon = 0.0000001;
  Parameters currentRates = startingParameters;
  updateLL(currentRates, evaluations);
  unsigned int llComputationsGrad = 0;
  unsigned int llComputationsLine = 0;
  unsigned int dimensions = startingParameters.dimensions();
  Parameters gradient(dimensions);
  do {
    std::vector<Parameters> closeRates(dimensions, currentRates);
    for (unsigned int i = 0; i < dimensions; ++i) {
      Parameters closeRates = currentRates;
      closeRates[i] += epsilon;
      updateLL(closeRates, evaluations);
      llComputationsGrad++;
      gradient[i] = (currentRates.getScore() - closeRates.getScore()) / (-epsilon);
    }
  } while (lineSearchParameters(evaluations, currentRates, gradient, llComputationsLine));
  return currentRates;
}

Parameters DTLOptimizer::optimizeParametersGlobalDTL(PerCoreGeneTrees &geneTrees, 
    PLLRootedTree &speciesTree, 
    RecModel model,
    Parameters *startingParameters)
{
  std::vector<Parameters> startingRates;
  if (startingParameters) {
    startingRates.push_back(*startingParameters);
  }
  else if (Enums::freeParameters(model) == 2) {
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
  Logger::info << "best: " << best << std::endl;
  return best; 
}


Parameters DTLOptimizer::optimizeParametersPerSpecies(PerCoreGeneTrees &geneTrees, 
    PLLRootedTree &speciesTree, 
    RecModel model)
{
  Parameters globalRates = optimizeParametersGlobalDTL(geneTrees, speciesTree, model);
  Parameters startingSpeciesRates(speciesTree.getNodesNumber(), globalRates);
  Parameters rates = DTLOptimizer::optimizeParameters(geneTrees, speciesTree, model, startingSpeciesRates);
  return rates; 
}




void DTLOptimizer::buildEvaluations(PerCoreGeneTrees &geneTrees, PLLRootedTree &speciesTree, RecModel recModel, Evaluations &evaluations)
{
  auto &trees = geneTrees.getTrees();
  evaluations.resize(trees.size());
  for (unsigned int i = 0; i < trees.size(); ++i) {
    auto &tree = trees[i];
    evaluations[i] = std::make_shared<ReconciliationEvaluation>(speciesTree, *tree.geneTree, tree.mapping, recModel, false);
  }
}


