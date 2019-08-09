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
    rates.ll = -std::numeric_limits<double>::infinity();
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
    ll = -std::numeric_limits<double>::infinity();
  }
  rates.setLL(ll);
}

void optimizeDTLRatesNewton(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, DTLRatesVector &ratesVector)
{
  std::vector<std::shared_ptr<ReconciliationEvaluation> > evaluations;
  //SubtreeRepeatsCache cache;
  for (auto &tree: geneTrees.getTrees()) {
    auto evaluation = std::make_shared<ReconciliationEvaluation> (speciesTree, tree.mapping, model, true);
    evaluations.push_back(evaluation);
   // cache.addTree(tree.tree, tree.mapping);
   // evaluation->getReconciliationModel()->setCache(&cache); 
  }
  double epsilon = -0.000001;
  unsigned int dimensions = Enums::freeParameters(model);
  unsigned int species = ratesVector.size();
  DTLRatesVector currentRates = ratesVector;
  updateLL(currentRates, geneTrees, evaluations);
  bool stop = false;
  while (!stop) {
    stop = true;
    DTLRatesVector gradient(species);
    for (unsigned int j = 0; j < species; ++j) {
      for (unsigned int i = 0; i < dimensions; ++i) {
        DTLRatesVector closeRates = currentRates;
        closeRates.getRates(j).rates[i] -= epsilon;
        updateLL(closeRates, geneTrees, evaluations);
        gradient.getRates(j).rates[i] = (currentRates.getLL() - closeRates.getLL()) / epsilon;
      }
    }
    double alpha = 0.1;
    while (alpha > 0.001) {
      gradient.normalize(alpha);
      DTLRatesVector proposal = currentRates + (gradient * alpha);
      updateLL(proposal, geneTrees, evaluations);
      if (currentRates.getLL() < proposal.getLL()) {
        currentRates = proposal;
        alpha *= 1.25;
        stop = false;
      } else {
        alpha *= 0.75;
      }
    }
  }
  Logger::info << "Newton ll: " << currentRates.getLL() << std::endl;
  ratesVector = currentRates;
}

DTLRatesVector DTLOptimizer::optimizeDTLRatesVector(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, DTLRatesVector *previousVector)
{
  unsigned int speciesNumber = speciesTree->inner_count + speciesTree->tip_count;
  DTLRatesVector starting(speciesNumber, DTLRates(1.0, 1.0, 1.0));
  if (previousVector && previousVector->size() == speciesNumber) {
    starting = *previousVector;    
  }
  optimizeDTLRatesNewton(geneTrees, speciesTree, model, starting);
  return starting;
  std::vector<std::shared_ptr<ReconciliationEvaluation> > evaluations;
  //SubtreeRepeatsCache cache;
  for (auto &tree: geneTrees.getTrees()) {
    auto evaluation = std::make_shared<ReconciliationEvaluation> (speciesTree, tree.mapping, model, true);
    evaluations.push_back(evaluation);
   // cache.addTree(tree.tree, tree.mapping);
   // evaluation->getReconciliationModel()->setCache(&cache); 
  }
  std::default_random_engine generator;
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
    //Logger::info << "ll " << rates[0].getLL() << std::endl;
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


DTLRates optimizeDTLRatesAux(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, 
    double dupRadius, double lossRadius, double transferRadius)
{
  std::vector<DTLRates> rates;
  if (Enums::accountsForTransfers(model)) {
    rates.push_back(DTLRates(dupRadius / 100.0, lossRadius / 100.0, 0.0));
    rates.push_back(DTLRates(dupRadius        , lossRadius        , 0.0));
    rates.push_back(DTLRates(dupRadius / 100.0, lossRadius        , transferRadius));
    rates.push_back(DTLRates(dupRadius        , lossRadius / 100.0, transferRadius));
  } else {
    rates.push_back(DTLRates(dupRadius / 100.0, lossRadius / 100, 0.0));
    rates.push_back(DTLRates(dupRadius        , lossRadius / 100, 0.0));
    rates.push_back(DTLRates(dupRadius / 100.0, lossRadius      , 0.0));
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
    double stepSize = 1 / double(iterations - 1);
    for (int i = 0; i < iterations; i++) {
      previous = current;
      current = x1 + ((x2 - x1) * i * stepSize);
      updateLL(current, geneTrees, speciesTree, model);
      llCalls++;
      if (current < bestRates) {
        bestRates = current;
      }
    } //while (downgrades < 2);
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

DTLRates optimizeDTLRatesNewtoon(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model, const DTLRates &startingRates)
{
  double epsilon = -0.000001;
  unsigned int dimensions = Enums::freeParameters(model);
  DTLRates currentRates;
  for (unsigned int i = 0; i < dimensions; ++i) {
    currentRates.rates[i] = startingRates.rates[i];
  }
  updateLL(currentRates, geneTrees, speciesTree, model);
  bool stop = false;
  unsigned int llComputations = 0;
  double minImprovement = 0.1;
  double lastAlpha = 0.1;
  while (!stop) {
    stop = true;
    DTLRates gradient;
    std::vector<DTLRates> closeRates(dimensions, currentRates);
    for (unsigned int i = 0; i < dimensions; ++i) {
      DTLRates closeRates = currentRates;
      closeRates.rates[i] -= epsilon;
      updateLL(closeRates, geneTrees, speciesTree, model);
      llComputations++;
      gradient.rates[i] = (currentRates.ll - closeRates.ll) / epsilon;
    }
    double alpha = lastAlpha / 8.0;
    lastAlpha = 0.0;
    while (alpha > 0.01) {
      gradient.normalize(alpha);
      DTLRates proposal = currentRates + (gradient * alpha);
      updateLL(proposal, geneTrees, speciesTree, model);
      llComputations++;
      if (currentRates.ll + minImprovement < proposal.ll) {
        currentRates = proposal;
        lastAlpha += alpha;
        alpha *= 2.0;
        stop = false;
      } else {
        alpha *= 0.5;
      }
    }
  }
  Logger::info << "Newton rates: " << currentRates << "(after " << llComputations << " recomputations) " << std::endl;
  return currentRates;
}

DTLRates DTLOptimizer::optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model)
{
  std::vector<DTLRates> startingRates;
  startingRates.push_back(DTLRates(1.0, 1.0, 1.0));
  startingRates.push_back(DTLRates(0.1, 0.2, 0.0));
  startingRates.push_back(DTLRates(0.02, 0.01, 0.01));
  startingRates.push_back(DTLRates(0.3, 0.3, 0.0));
  DTLRates best;
  best.ll = -10000000000;
  for (auto rates: startingRates) {
    DTLRates newRates = optimizeDTLRatesNewtoon(geneTrees, speciesTree, model, rates);
    bool stop = (fabs(newRates.ll - best.ll) < 1.0);
    if (newRates.ll > best.ll) {
      best = newRates;
    }
    if (stop) {
      break;
    }
  }
  return best; 
  std::vector<DTLRates> bestRates;
  if (Enums::accountsForTransfers(model)) {
    bestRates.push_back(optimizeDTLRatesAux(geneTrees, speciesTree, model, 1.0, 1.0, 1.0));
    bestRates.push_back(optimizeDTLRatesAux(geneTrees, speciesTree, model, 0.2, 0.4, 0.2));
  } else {
    bestRates.push_back(optimizeDTLRatesAux(geneTrees, speciesTree, model, 1.0, 1.0, 0.0));
    bestRates.push_back(optimizeDTLRatesAux(geneTrees, speciesTree, model, 0.2, 0.2, 0.0));
  }
  sort(bestRates.begin(), bestRates.end());

  return bestRates[0];
}
  
