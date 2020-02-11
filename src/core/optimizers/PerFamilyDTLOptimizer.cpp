#include <optimizers/PerFamilyDTLOptimizer.hpp>

#include <parallelization/ParallelContext.hpp>
#include <trees/JointTree.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <IO/Logger.hpp>
#include <limits>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <optimizers/DTLOptimizer.hpp>

static bool isValidLikelihood(double ll) {
  return std::isnormal(ll) && ll < -0.0000001;
}

static void updateLL(Parameters &rates, JointTree &jointTree) {
  rates.ensurePositivity();
  jointTree.setRates(rates);
  rates.setScore(jointTree.computeReconciliationLoglk());
  if (!isValidLikelihood(rates.getScore())) {
    rates.setScore(-100000000000.0);
  }
}

void PerFamilyDTLOptimizer::findBestRatesDL(JointTree &jointTree,
    double minDup, double maxDup,
    double minLoss, double maxLoss, unsigned int steps,
    double &bestDup,
    double &bestLoss,
    double &bestLL) 
{
  bestLL = std::numeric_limits<double>::lowest();
  auto totalSteps = steps * steps;
  auto begin = ParallelContext::getBegin(totalSteps);
  auto end = ParallelContext::getEnd(totalSteps);
  for (auto s = begin; s < end; ++s) {
    auto i = s / steps;
    auto j = s % steps;
    double dup = minDup + (maxDup - minDup) * double(i) / double(steps);
    double loss = minLoss + (maxLoss - minLoss) * double(j) / double(steps);
    Parameters rates(dup, loss);
    rates.ensurePositivity();
    double newLL = jointTree.computeReconciliationLoglk();
    if (!isValidLikelihood(newLL)) {
      continue;
    }
    if (newLL > bestLL) { 
      bestDup = rates[0];
      bestLoss = rates[1];
      bestLL = newLL;
    }
  }
  unsigned int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  jointTree.setRates(Parameters(bestDup, bestLoss));
}

void PerFamilyDTLOptimizer::findBestRatesDTL(JointTree &jointTree,
    double minDup, double maxDup,
    double minLoss, double maxLoss, 
    double minTrans, double maxTrans, 
    unsigned int steps,
    double &bestDup,
    double &bestLoss,
    double &bestTrans,
    double &bestLL) 
{
  bestLL = std::numeric_limits<double>::lowest();
  unsigned int totalSteps = steps * steps * steps;
  auto begin = ParallelContext::getBegin(totalSteps);
  auto end = ParallelContext::getEnd(totalSteps);
  for (auto s = begin; s < end; ++s) {
    auto i = s / (steps * steps);
    auto j = (s / steps) % steps;
    auto k = s % steps;
    double dup = minDup + (maxDup - minDup) * double(i) / double(steps);
    double loss = minLoss + (maxLoss - minLoss) * double(j) / double(steps);
    double trans = minTrans + (maxTrans - minTrans) * double(k) / double(steps);
    Parameters rates(dup, loss, trans);
    rates.ensurePositivity();
    jointTree.setRates(rates);
    double newLL = jointTree.computeReconciliationLoglk();
    if (!isValidLikelihood(newLL)) {
      continue;
    }
    if (newLL > bestLL) { 
      bestDup = rates[0];
      bestLoss = rates[1];
      bestTrans = rates[2];
      bestLL = newLL;
    }
  }
  unsigned int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  ParallelContext::broadcastDouble(bestRank, bestTrans);
  jointTree.setRates(Parameters(bestDup, bestLoss, bestTrans));
}


void PerFamilyDTLOptimizer::optimizeDTLRates(JointTree &jointTree, RecOpt method) {
  switch(method) {
  case RecOpt::Grid:
    optimizeDTLRatesWindow(jointTree);
    break;
  case RecOpt::Simplex:
    optimizeRateSimplex(jointTree, true);
    break;
  case RecOpt::Gradient:
    optimizeDTLRatesGradient(jointTree);
    break;
  }
}

void PerFamilyDTLOptimizer::optimizeDLRates(JointTree &jointTree, RecOpt method) {
  switch(method) {
  case RecOpt::Grid:
    optimizeDLRatesWindow(jointTree);
    break;
  case RecOpt::Simplex:
    optimizeRateSimplex(jointTree, false);
    break;
  case RecOpt::Gradient:
    optimizeDTLRatesGradient(jointTree);
    break;
  }
}


void PerFamilyDTLOptimizer::optimizeDLRatesWindow(JointTree &jointTree) {
  Logger::timed << "Start optimizing DL rates" << std::endl;
  
  double bestLL = std::numeric_limits<double>::lowest();
  double newLL = 0;
  double bestDup = 0.0;
  double bestLoss = 0.0;
  double minDup = 0.0;
  double maxDup = 10.0;
  double minLoss = 0.0;
  double maxLoss = 10.0;
  unsigned int steps = 10;
  double epsilon = 0.001;
  do {
    bestLL = newLL;
    findBestRatesDL(jointTree, minDup, maxDup, minLoss, maxLoss, steps, bestDup, bestLoss, newLL);
    double offsetDup = (maxDup - minDup) / steps;
    double offsetLoss =(maxLoss - minLoss) / steps;
    minDup = std::max(0.0, bestDup - offsetDup);
    maxDup = bestDup + offsetDup;
    minLoss = std::max(0.0, bestLoss - offsetLoss);
    maxLoss = bestLoss + offsetLoss;
  } while (fabs(newLL - bestLL) > epsilon);
  Logger::info << " best rates: " << bestDup << " " << bestLoss <<  " " << newLL << std::endl;
  if  (!isValidLikelihood(newLL)) {
    Logger::error << "Invalid likelihood " << newLL << std::endl;
    ParallelContext::abort(10);
  }
}
    
void PerFamilyDTLOptimizer::optimizeDTLRatesWindow(JointTree &jointTree) {
  Logger::timed << "Start optimizing DTL rates" << std::endl;
  double bestLL = std::numeric_limits<double>::lowest();
  double newLL = 0;
  double bestDup = 0.0;
  double bestLoss = 0.0;
  double bestTrans = 0.0;
  double minDup = 0.0;
  double maxDup = 1.0;
  double minLoss = 0.0;
  double maxLoss = 1.0;
  double minTrans = 0.0;
  double maxTrans = 1.0;
  unsigned int steps = 5;
  double epsilon = 0.01;
  do {
    bestLL = newLL;
    findBestRatesDTL(jointTree, minDup, maxDup, minLoss, maxLoss, minTrans, maxTrans, steps, bestDup, bestLoss, bestTrans, newLL);
    double offsetDup = (maxDup - minDup) / steps;
    double offsetLoss = (maxLoss - minLoss) / steps;
    double offsetTrans = (maxTrans - minTrans) / steps;
    minDup = std::max(0.0, bestDup - offsetDup);
    maxDup = bestDup + offsetDup;
    minLoss = std::max(0.0, bestLoss - offsetLoss);
    maxLoss = bestLoss + offsetLoss;
    minTrans = std::max(0.0, bestTrans - offsetTrans);
    maxTrans = bestTrans + offsetTrans;
  } while (fabs(newLL - bestLL) > epsilon);
  if  (!isValidLikelihood(newLL)) {
    Logger::error << "Invalid likelihood " << newLL << std::endl;
    ParallelContext::abort(10);
  }
}



/*
 *  Find the point between  r1 and r2. Parallelized over the iterations
 */
Parameters findBestPoint(Parameters r1, Parameters r2, unsigned int iterations, JointTree &jointTree) 
{
  Parameters best = r1;
  best.setScore(-100000000000);
  unsigned int bestI = 0;
  for (auto i = ParallelContext::getRank(); i < iterations; i += ParallelContext::getSize()) {
    Parameters current = r1 + ((r2 - r1) * (double(i) / double(iterations - 1)));
    updateLL(current, jointTree);
    if (current < best) {
      best = current;
      bestI = i;
    }
  }
  unsigned int bestRank = 0;
  double ll = best.getScore();
  ParallelContext::getMax(ll, bestRank);
  ParallelContext::broadcastUInt(bestRank, bestI);
  if (ParallelContext::getRank() != bestRank) {
    best = r1 + ((r2 - r1) * (double(bestI) / double(iterations - 1)));
    best.ensurePositivity();
  }
  ParallelContext::broadcastDouble(bestRank, ll);
  best.setScore(ll);
  return best;
}



void PerFamilyDTLOptimizer::optimizeRateSimplex(JointTree &jointTree, bool transfers)
{
  Logger::timed << "Starting DTL rates optimization" << std::endl;
  std::vector<Parameters> rates;
  rates.push_back(Parameters(0.01, 0.01, 0.01));
  rates.push_back(Parameters(1.0, 0.01, 0.01));
  rates.push_back(Parameters(0.01, 1.0, 1.0));
  if (transfers) {
    rates.push_back(Parameters(0.01, 0.01, 1.0));
  }
  for (auto &r: rates) {
    updateLL(r, jointTree);
  }
  Parameters worstRate;
  unsigned int currentIt = 0;
  while (worstRate.distance(rates.back()) > 0.005) {
    sort(rates.begin(), rates.end());
    worstRate = rates.back();
    // centroid
    Parameters x0;
    for (unsigned int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion, exansion and contraction at the same time
    Parameters x1 = x0 - (x0 - rates.back()) * 0.5;  
    Parameters x2 = x0 + (x0 - rates.back()) * 1.5;  
    unsigned int iterations = 8;
    Parameters xr = findBestPoint(x1, x2, iterations, jointTree);
    if (xr < rates[rates.size() - 1] ) {
      rates.back() = xr;
    }
    currentIt++;
  }
  sort(rates.begin(), rates.end());
  updateLL(rates[0], jointTree);
  Logger::timed << "Simplex converged after " << currentIt << " iterations" << std::endl;
  Logger::timed << rates[0] << std::endl;
}

  
void PerFamilyDTLOptimizer::optimizeDTLRatesGradient(JointTree &jointTree)
{
  Evaluations evaluations;
  evaluations.push_back(jointTree.getReconciliationEvaluationPtr());
  Parameters rates = DTLOptimizer::optimizeParameters(evaluations, jointTree.getRatesVector());
  Logger::info << "Per family rates: " << rates << std::endl;
  jointTree.setRates(rates);
}

