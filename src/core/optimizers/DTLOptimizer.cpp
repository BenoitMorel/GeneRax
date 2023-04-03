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
  Logger::info << "lineSearch " << currentRates.getScore() <<  " gradient: " << gradient << std::endl;
  //Logger::info << "minimprov " << settings.lineSearchMinImprovement <<  std::endl; 
  while (alpha > minAlpha) {
    currentGradient.normalize(alpha);
    Parameters proposal = currentRates + (currentGradient * alpha);
    function.evaluate(proposal);
    llComputationsLine++;
    if (currentRates.getScore() + settings.lineSearchMinImprovement
        < proposal.getScore()) {
      Logger::info << "Improv alpha=" << alpha << " score=" << proposal.getScore() << " p=" << proposal << std::endl;
      currentRates = proposal;
      noImprovement = false;
      alpha *= 1.5;
    } else {
      alpha *= 0.5;
      //Logger::info << std::setprecision(15) << "No improv alpha=" << alpha << " score_to_beat=" << currentRates.getScore() << " p=" << proposal  << std::endl;
      if (!noImprovement) {
        return true;
      }
    }
  }
  return !noImprovement;
}


Parameters optimizeParametersGradient(FunctionToOptimize &function,
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
  Logger::info << "Computing gradient..." << std::endl;
  do {
    std::vector<Parameters> closeRates(dimensions, currentRates);
    for (unsigned int i = 0; i < dimensions; ++i) {
      Parameters closeRates = currentRates;
      closeRates[i] += epsilon;
      //Logger::info << "Close rates: " << closeRates << std::endl;
      function.evaluate(closeRates);
      llComputationsGrad++;
      gradient[i] = (currentRates.getScore() - closeRates.getScore()) / (-epsilon);
      //Logger::info << " GRAD: currll=" << currentRates.getScore() << " newll=" << closeRates.getScore() << " diff=" << currentRates.getScore() - closeRates.getScore() << " grad[i]=" << gradient[i] << " epsilon=" << epsilon  << std::endl;
    }
  } while (lineSearchParameters(function, currentRates, gradient, llComputationsLine, settings));
  function.evaluate(currentRates);
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


static Parameters findBestPointNelderMear(Parameters r1, 
    Parameters r2, 
    unsigned int iterations, 
    FunctionToOptimize &function) 
{
  Parameters best = r1;
  best.setScore(-100000000000);
  unsigned int bestI = 0;
  for (auto i = 0; i < iterations; ++i) {
    Parameters current = r1 + ((r2 - r1) * (double(i) / double(iterations - 1)));
    function.evaluate(current);
    if (current < best) {
      best = current;
      bestI = i;
    }
  }
  return best;
}



static Parameters optimizeParametersNelderMear(FunctionToOptimize &function, 
    const Parameters &startingParameters,
    OptimizationSettings settings = OptimizationSettings())
{
  std::vector<Parameters> rates;
  rates.push_back(startingParameters);
  auto N = startingParameters.dimensions();
  for (unsigned int r = 0; r < N; ++r) {
    auto p = startingParameters;
    p[r] -=  0.09;
    rates.push_back(p);
  }
  // n + 1 points for n dimensions
  assert(rates.size() == N + 1);

  for (auto &r: rates) {
    function.evaluate(r);
  }
  Parameters worstRate = startingParameters;
  unsigned int currentIt = 0;
 
  while (worstRate.distance(rates.back()) > 0.005) {
    std::sort(rates.begin(), rates.end());
    worstRate = rates.back();
    // centroid
    Parameters x0(worstRate.dimensions());
    for (unsigned int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion, exansion and contraction at the same time
    Parameters x1 = x0 - (x0 - rates.back()) * 0.5;  
    Parameters x2 = x0 + (x0 - rates.back()) * 1.5;  
    unsigned int iterations = 8;
    Parameters xr = findBestPointNelderMear(x1, x2, iterations, function);
    if (xr < rates[rates.size() - 1] ) {
      rates.back() = xr;
    }
    currentIt++;
  }
  std::sort(rates.begin(), rates.end());
  function.evaluate(rates[0]);
  Logger::timed << "Simplex converged after " << currentIt << " iterations" << std::endl;
  return rates[0];
}

Parameters DTLOptimizer::optimizeParameters(FunctionToOptimize &function,
    Parameters startingParameters,
    OptimizationSettings settings)
{
  switch(settings.strategy) {
  case RecOpt::Gradient:
      return optimizeParametersGradient(function, 
          startingParameters, 
          settings);
  case RecOpt::Simplex:
      return optimizeParametersNelderMear(function, 
          startingParameters,
          settings);
  default:
      assert(false);
      return startingParameters;
  }
}
      

static Parameters optimizeParametersCorax(FunctionToOptimize &function, 
    const Parameters &startingParameters,
    OptimizationSettings settings = OptimizationSettings())
{
  unsigned int xnum = startingParameters.dimensions();
  //double *xparameters = &startingParameters.getVector()[0];
  //double **x = &xparameters;
  //Parameters paramMin(

  return startingParameters;
}

/*
CORAX_EXPORT double corax_opt_minimize_lbfgsb_multi(
    unsigned int  xnum,
    double      **x,
    double      **xmin,
    double      **xmax,
    int         **bound,
    unsigned int *n,
    unsigned int  nmax,
    double        factr,
    double        pgtol,
    void         *params,
    double (*target_funk)(void *, double **, double *, int *));
*/
