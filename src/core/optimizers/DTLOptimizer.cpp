#include <optimizers/DTLOptimizer.hpp>
#include <iomanip>
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <limits>
#include <algorithm>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <iostream>
#include <cmath>

static bool isValidLikelihood(double ll) {
  return std::isnormal(ll) && ll < -0.0000001;
}



static bool lineSearchParameters(FunctionToOptimize &function,
    Parameters &currentRates, 
    const Parameters &gradient, 
    unsigned int &llComputationsLine,
    const OptimizationSettings &settings
    )
{
  double alpha = 0.1; 
  const double minAlpha = settings.minAlpha;
  Parameters currentGradient(gradient);
  bool noImprovement = true;
  //Logger::info << "lineSearch " << currentRates.getScore() <<  " gradient: " << gradient << std::endl;
  while (alpha > minAlpha) {
    currentGradient.normalize(alpha);
    Parameters proposal = currentRates + (currentGradient * alpha);
    function.evaluate(proposal);
    llComputationsLine++;
    if (currentRates.getScore() + settings.lineSearchMinImprovement
        < proposal.getScore()) {
      //Logger::info << "Improv alpha=" << alpha << " score=" << proposal.getScore() << " p=" << proposal << std::endl;
      currentRates = proposal;
      noImprovement = false;
      alpha *= 1.5;
    } else {
      alpha *= 0.5;
      if (!noImprovement) {
        return true;
      }
      //Logger::info << std::setprecision(15) << "No improv alpha=" << alpha << " score=" << proposal.getScore() << " p=" << proposal  << std::endl;
    }
  }
  return !noImprovement;
}


Parameters DTLOptimizer::optimizeParameters(FunctionToOptimize &function,
    const Parameters &startingParameters,
    OptimizationSettings settings)
{
  if (startingParameters.dimensions() == 0) {
    return Parameters();
  }
  double epsilon = settings.epsilon;
  Parameters currentRates = startingParameters;
  function.evaluate(currentRates);
  unsigned int llComputationsGrad = 0;
  unsigned int llComputationsLine = 0;
  unsigned int dimensions = startingParameters.dimensions();
  Parameters gradient(dimensions);
  do {
    std::vector<Parameters> closeRates(dimensions, currentRates);
    for (unsigned int i = 0; i < dimensions; ++i) {
      Parameters closeRates = currentRates;
      closeRates[i] += epsilon;
      function.evaluate(closeRates);
      llComputationsGrad++;
      gradient[i] = (currentRates.getScore() - closeRates.getScore()) / (-epsilon);
      //Logger::info << "gradient" << currentRates.getScore() - closeRates.getScore() << " " << gradient[i] << std::endl;
    }
  } while (lineSearchParameters(function, currentRates, gradient, llComputationsLine, settings));
  /*
  Parameters minusGradient(gradient);
  for (unsigned int i = 0; i < gradient.dimensions(); ++i) {
    minusGradient[i] = -gradient[i];
  }
  Logger::info << "Final line search with minus gradient" << std::endl;
  lineSearchParameters(function, currentRates, minusGradient, llComputationsLine, settings);
  */
  return currentRates;
}

class PerCoreFunction: public FunctionToOptimize {
public:
  PerCoreFunction(PerCoreEvaluations &evaluations):_evaluations(evaluations){}
  virtual double evaluate(Parameters &parameters) {
    parameters.ensurePositivity();
    double ll = 0.0;
    for (auto evaluation: _evaluations) {
      evaluation->setRates(parameters);
      ll += evaluation->evaluate();
    }
    ParallelContext::sumDouble(ll);
    if (!isValidLikelihood(ll)) {
      ll = -std::numeric_limits<double>::infinity();
    }
    parameters.setScore(ll);
    return ll;
  }
private:
  PerCoreEvaluations &_evaluations;
};

Parameters DTLOptimizer::optimizeParameters(PerCoreEvaluations &evaluations,
    const Parameters &startingParameters,
    OptimizationSettings settings)
{
  PerCoreFunction function(evaluations);
  return optimizeParameters(function, startingParameters, settings);
}

ModelParameters DTLOptimizer::optimizeModelParameters(PerCoreEvaluations &evaluations,
    bool optimizeFromStartingParameters,
    const ModelParameters &startingParameters,
    OptimizationSettings settings)
{

  ModelParameters res = startingParameters;
  if (!startingParameters.info.perFamilyRates) {
    const Parameters *startingRates = optimizeFromStartingParameters ? &startingParameters.rates :  nullptr;
    res.rates = DTLOptimizer::optimizeParametersGlobalDTL(evaluations, startingRates, settings);
  } else {
    ParallelContext::pushSequentialContext(); // work locally
    for (unsigned int i = 0; i < evaluations.size(); ++i) {
      Parameters localRates = startingParameters.getRates(i);
      const Parameters *startingRates = optimizeFromStartingParameters ? &localRates : nullptr;
      PerCoreEvaluations localEvaluation;
      localEvaluation.push_back(evaluations[i]);
      localRates = DTLOptimizer::optimizeParametersGlobalDTL(localEvaluation, startingRates, settings);
      res.setRates(i, localRates);
    }
    ParallelContext::popContext();
  }
  return res;
}


Parameters DTLOptimizer::optimizeParametersGlobalDTL(PerCoreEvaluations &evaluations, 
    const Parameters *startingParameters,
    OptimizationSettings settings)
{
  unsigned int freeParameters = 0;
  if (evaluations.size()) {
    freeParameters = Enums::freeParameters(evaluations[0]->getRecModel());
  }
  ParallelContext::maxUInt(freeParameters);
  if (freeParameters == 0) {
    return Parameters();
  }
  std::vector<Parameters> startingRates;
  if (startingParameters) {
    startingRates.push_back(*startingParameters);
  }
  if (freeParameters == 1) {
    Parameters p(1);
    p[0] = 0.1;
    startingRates.push_back(p);
    p[0] = 0.3;
    startingRates.push_back(p);
    p[0] = 1.0;
    startingRates.push_back(p);
    p[0] = 10.0;
    startingRates.push_back(p);
  } else if (freeParameters == 2) {
    startingRates.push_back(Parameters(0.1, 0.2));
    startingRates.push_back(Parameters(0.2, 0.2));
    startingRates.push_back(Parameters(0.5, 0.5));
    startingRates.push_back(Parameters(0.5, 1.0));
    startingRates.push_back(Parameters(0.01, 0.01));
  } else if (freeParameters == 3) {
    startingRates.push_back(Parameters(0.1, 0.2, 0.1));
    startingRates.push_back(Parameters(0.01, 0.01, 0.01));
  } else {
    startingRates.push_back(Parameters(0.5, 0.5, 0.2, 0.0));
    startingRates.push_back(Parameters(0.1, 0.2, 0.1, 0.1));
    startingRates.push_back(Parameters(0.2, 0.2, 0.0, 0.1));
    startingRates.push_back(Parameters(0.01, 0.01, 0.01, 0.01));
  }
  ParallelContext::barrier();
  Parameters best;
  best.setScore(-10000000000);
  for (auto rates: startingRates) {
    Parameters newRates = optimizeParameters(evaluations, rates, settings);
    bool stop = (fabs(newRates.getScore() - best.getScore()) 
        < settings.optimizationMinImprovement);
    stop = false;
    if (newRates.getScore() > best.getScore()) {
      best = newRates;
    }
    if (stop) {
      break;
    }
  }
  return best; 
}



Parameters DTLOptimizer::optimizeParametersPerSpecies(PerCoreEvaluations &evaluations, unsigned int speciesNodesNumber) 
{
  Parameters globalRates = optimizeParametersGlobalDTL(evaluations);
  Parameters startingSpeciesRates(speciesNodesNumber, globalRates);
  Parameters rates = DTLOptimizer::optimizeParameters(evaluations, startingSpeciesRates);
  return rates; 
}






