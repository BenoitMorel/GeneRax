#include <optimizers/DTLOptimizer.hpp>

#include <likelihoods/SubtreeRepeatsCache.hpp>
#include <ParallelContext.hpp>
#include <trees/JointTree.hpp>
#include <IO/Logger.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <limits>
#include <algorithm>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <iostream>
#include <cmath>

bool isValidLikelihood(double ll) {
  return std::isnormal(ll) && ll < -0.0000001;
}

void updateLL(DTLRates &rates, JointTree &jointTree) {
  rates.ensureValidity();
  jointTree.setRates(DTLRatesVector(rates));
  rates.ll = jointTree.computeReconciliationLoglk();
  if (!isValidLikelihood(rates.ll)) {
    rates.ll = -100000000000.0;
  }
}

void updateLL(DTLRates &rates, PerCoreGeneTrees &trees, pll_rtree_t *speciesTree, RecModel model) {
  rates.ensureValidity();
  rates.ll = 0.0;
  for (auto &tree: trees.getTrees()) {
    ReconciliationEvaluation evaluation(speciesTree, tree.mapping, model, true);
    evaluation.setRates(rates.rates[0], rates.rates[1], rates.rates[2]);
    rates.ll += evaluation.evaluate(tree.tree);
  }
  ParallelContext::sumDouble(rates.ll);
  if (!isValidLikelihood(rates.ll)) {
    rates.ll = -100000000000.0;
  }
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
    ll = -100000000000.0;
  }
  rates.setLL(ll);
}


void DTLOptimizer::findBestRatesDL(JointTree &jointTree,
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
    jointTree.setRates(DTLRatesVector(DTLRates(dup, loss)));
    double newLL = jointTree.computeReconciliationLoglk();
    if (!isValidLikelihood(newLL)) {
      continue;
    }
    if (newLL > bestLL) { 
      bestDup = dup;
      bestLoss = loss;
      bestLL = newLL;
    }
  }
  unsigned int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  jointTree.setRates(DTLRatesVector(DTLRates(bestDup, bestLoss)));
}

void DTLOptimizer::findBestRatesDTL(JointTree &jointTree,
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
    jointTree.setRates(DTLRatesVector(DTLRates(dup, loss, trans)));
    double newLL = jointTree.computeReconciliationLoglk();
    if (!isValidLikelihood(newLL)) {
      continue;
    }
    if (newLL > bestLL) { 
      bestDup = dup;
      bestLoss = loss;
      bestTrans = trans;
      bestLL = newLL;
    }
  }
  unsigned int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  ParallelContext::broadcastDouble(bestRank, bestTrans);
  jointTree.setRates(DTLRatesVector(DTLRates(bestDup, bestLoss, bestTrans)));
}


void DTLOptimizer::optimizeDTLRates(JointTree &jointTree, RecOpt method) {
  switch(method) {
  case Grid:
    optimizeDTLRatesWindow(jointTree);
    break;
  case Simplex:
    optimizeRateSimplex(jointTree, true);
    break;
  }
  /*
  Logger::info << "best rates " << std::endl;
  Logger::info << "D " << jointTree.getDupRate() << std::endl;
  Logger::info << "L " << jointTree.getLossRate() << std::endl;
  Logger::info << "T " << jointTree.getTransferRate() << std::endl;
  */
}

void DTLOptimizer::optimizeDLRates(JointTree &jointTree, RecOpt method) {
  switch(method) {
  case Grid:
    optimizeDLRatesWindow(jointTree);
    break;
  case Simplex:
    optimizeRateSimplex(jointTree, false);
    break;
  }
  /*
  Logger::info << "best rates " << std::endl;
  Logger::info << "D " << jointTree.getDupRate() << std::endl;
  Logger::info << "L " << jointTree.getLossRate() << std::endl;
  */
}


void DTLOptimizer::optimizeDLRatesWindow(JointTree &jointTree) {
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
    
void DTLOptimizer::optimizeDTLRatesWindow(JointTree &jointTree) {
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
DTLRates findBestPoint(DTLRates r1, DTLRates r2, unsigned int iterations, JointTree &jointTree) 
{
  DTLRates best = r1;
  best.ll = -100000000000;
  unsigned int bestI = 0;
  for (auto i = ParallelContext::getRank(); i < iterations; i += ParallelContext::getSize()) {
    DTLRates current = r1 + ((r2 - r1) * (double(i) / double(iterations - 1)));
    updateLL(current, jointTree);
    if (current < best) {
      best = current;
      bestI = i;
    }
  }
  unsigned int bestRank = 0;
  ParallelContext::getMax(best.ll, bestRank);
  ParallelContext::broadcastUInt(bestRank, bestI);
  if (ParallelContext::getRank() != bestRank) {
    best = r1 + ((r2 - r1) * (double(bestI) / double(iterations - 1)));
    best.ensureValidity();
  }
  ParallelContext::broadcastDouble(bestRank, best.ll);
  return best;
}



void DTLOptimizer::optimizeRateSimplex(JointTree &jointTree, bool transfers)
{
  Logger::timed << "Starting DTL rates optimization" << std::endl;
  std::vector<DTLRates> rates;
  rates.push_back(DTLRates(0.01, 0.01, 0.01));
  rates.push_back(DTLRates(1.0, 0.01, 0.01));
  rates.push_back(DTLRates(0.01, 1.0, 1.0));
  if (transfers) {
    rates.push_back(DTLRates(0.01, 0.01, 1.0));
  }
  for (auto &r: rates) {
    updateLL(r, jointTree);
  }
  DTLRates worstRate;
  unsigned int currentIt = 0;
  while (worstRate.distance(rates.back()) > 0.005) {
    sort(rates.begin(), rates.end());
    worstRate = rates.back();
    // centroid
    DTLRates x0;
    for (unsigned int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion, exansion and contraction at the same time
    DTLRates x1 = x0 - (x0 - rates.back()) * 0.5;  
    DTLRates x2 = x0 + (x0 - rates.back()) * 1.5;  
    unsigned int iterations = 8;
    DTLRates xr = findBestPoint(x1, x2, iterations, jointTree);
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

DTLRatesVector DTLOptimizer::optimizeDTLRatesVector(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, DTLRatesVector *previousVector)
{
  std::vector<std::shared_ptr<ReconciliationEvaluation> > evaluations;
  SubtreeRepeatsCache cache;
  for (auto &tree: geneTrees.getTrees()) {
    auto evaluation = std::make_shared<ReconciliationEvaluation> (speciesTree, tree.mapping, model, true);
    evaluations.push_back(evaluation);
    cache.addTree(tree.tree, tree.mapping);
    evaluation->getReconciliationModel()->setCache(&cache); 
  }
  std::default_random_engine generator;
  unsigned int speciesNumber = speciesTree->inner_count + speciesTree->tip_count;
  unsigned int freeRates = speciesNumber;
  if (Enums::accountsForTransfers(model)) {
    freeRates *= 3;
  } else {
    freeRates *= 2;
  }
  const double minRates = 0.00001;
  const double maxRates = 1.0;
  DTLRatesVector nullRates(speciesNumber);
  std::vector<DTLRatesVector> rates(freeRates + 1, nullRates);
  assert(rates.size() > 2);
  for (auto &rate: rates) {
    rate.initRandom(minRates, maxRates, generator);
  }
  rates[0] = nullRates;
  rates[1] = DTLRatesVector(speciesNumber, DTLRates(maxRates, maxRates, maxRates));
  if (previousVector && previousVector->size() == speciesNumber) {
    rates[2] = *previousVector;
  }
  int llCalls = 0;
  for (auto &r: rates) {
    updateLL(r, geneTrees, evaluations);
    llCalls++;
  }
  DTLRatesVector worstRate(nullRates);
  int currentIt = 0;
  while (worstRate.distance(rates.back()) > 0.001) {
    sort(rates.begin(), rates.end());
    Logger::info << "ll " << rates[0].getLL() << std::endl;
    worstRate = rates.back();
    // centroid
    DTLRatesVector x0(nullRates);
    for (unsigned int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion, exansion and contraction at the same time
    DTLRatesVector x1 = x0 - (x0 - rates.back()) * 0.5;  
    DTLRatesVector x2 = x0 + (x0 - rates.back()) * 1.5;  
    int iterations = 8;
    DTLRatesVector bestRates = x1;
    bestRates.setLL(-100000000000.0);
    DTLRatesVector previous(nullRates);
    DTLRatesVector current = x1;
    current.setLL(-10000000000.0);
    double stepSize = 1 / double(iterations - 1);
    for (int i = 0; i < iterations; i++) {
      previous = current;
      current = x1 + ((x2 - x1) * i * stepSize);
      updateLL(current, geneTrees, evaluations);
      llCalls++;
      if (current < bestRates) {
        bestRates = current;
      }
    }
    if (bestRates < rates[rates.size() - 1] ) {
      rates.back() = bestRates;
    }
    currentIt++;
  }
  Logger::info << "Best rates " << rates[0] << std::endl;
  ParallelContext::barrier();
  Logger::info << "end yoyo" << std::endl;
  sort(rates.begin(), rates.end());
  updateLL(rates[0], geneTrees, evaluations);
  llCalls++;
  return rates[0];
}

DTLRates DTLOptimizer::optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model)
{
  std::vector<DTLRates> rates;
  if (Enums::accountsForTransfers(model)) {
    rates.push_back(DTLRates(0.01, 0.01, 0.0));
    rates.push_back(DTLRates(1.0, 1.0, 0.0));
    rates.push_back(DTLRates(0.01, 1.0, 1.0));
    rates.push_back(DTLRates(0.5, 1.0, 0.5));
  } else {
    rates.push_back(DTLRates(0.01, 0.01, 0.0));
    rates.push_back(DTLRates(1.0, 0.01, 0.0));
    rates.push_back(DTLRates(0.01, 1.0, 0.0));
  }
  int llCalls = 0;
  for (auto &r: rates) {
    updateLL(r, geneTrees, speciesTree, model);
    llCalls++;
  }
  DTLRates worstRate;
  int currentIt = 0;
  while (worstRate.distance(rates.back()) > 0.005) {
    sort(rates.begin(), rates.end());
    worstRate = rates.back();
    // centroid
    DTLRates x0;
    for (unsigned int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion, exansion and contraction at the same time
    DTLRates x1 = x0 - (x0 - rates.back()) * 0.5;  
    DTLRates x2 = x0 + (x0 - rates.back()) * 1.5;  
    int iterations = 8;
    DTLRates bestRates = x1;
    bestRates.ll = -100000000000;
    DTLRates previous;
    DTLRates current = x1;
    current.ll = -10000000000;
    int i = 0;
    int downgrades = 0;
    int upgrades = 0;
    double stepSize = 1 / double(iterations - 1);
    do { //for (int i = 0; i < iterations; i++) {
      if (upgrades > 2) {
        stepSize *= 2;
      }
      previous = current;
      current = x1 + ((x2 - x1) * i * stepSize);
      updateLL(current, geneTrees, speciesTree, model);
      llCalls++;
      if (current < bestRates) {
        bestRates = current;
      }
      if (current < previous) {
        downgrades = 0;
        upgrades++;
      } else {
        upgrades = 0;
        downgrades ++;
      }
      ++i;
    } while (downgrades < 2);
    if (bestRates < rates[rates.size() - 1] ) {
      rates.back() = bestRates;
    }
    currentIt++;
  }
  sort(rates.begin(), rates.end());
  updateLL(rates[0], geneTrees, speciesTree, model);
  llCalls++;
  Logger::timed << "Simplexed converged to " << rates[0] << std::endl;
  //  Logger::timed << "Simplex converged after " << currentIt << " iterations and " << llCalls << " calls: " << rates[0] << std::endl;
  return rates[0];
}

  
