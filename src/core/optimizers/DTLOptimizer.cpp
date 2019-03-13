#include <optimizers/DTLOptimizer.hpp>

#include <ParallelContext.hpp>
#include <treeSearch/JointTree.hpp>
#include <IO/Logger.hpp>
#include <PerCoreGeneTrees.hpp>
#include <limits>
#include <algorithm>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <iostream>

bool isValidLikelihood(double ll) {
  return isnormal(ll) && ll < -0.0000001;
}

void DTLOptimizer::findBestRatesDL(JointTree &jointTree,
    double minDup, double maxDup,
    double minLoss, double maxLoss, int steps,
    double &bestDup,
    double &bestLoss,
    double &bestLL) 
{
  bestLL = numeric_limits<double>::lowest();
  int totalSteps = pow(steps, 2);
  int begin = ParallelContext::getBegin(totalSteps);
  int end = ParallelContext::getEnd(totalSteps);
  for (int s = begin; s < end; ++s) {
    int i = s / steps;
    int j = s % steps;
    double dup = minDup + (maxDup - minDup) * double(i) / double(steps);
    double loss = minLoss + (maxLoss - minLoss) * double(j) / double(steps);
    jointTree.setRates(dup, loss);
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
  int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  jointTree.setRates(bestDup, bestLoss);
}

void DTLOptimizer::findBestRatesDTL(JointTree &jointTree,
    double minDup, double maxDup,
    double minLoss, double maxLoss, 
    double minTrans, double maxTrans, 
    int steps,
    double &bestDup,
    double &bestLoss,
    double &bestTrans,
    double &bestLL) 
{
  bestLL = numeric_limits<double>::lowest();
  int totalSteps = pow(steps, 3);
  int begin = ParallelContext::getBegin(totalSteps);
  int end = ParallelContext::getEnd(totalSteps);
  for (int s = begin; s < end; ++s) {
    int i = s / (steps * steps);
    int j = (s / steps) % steps;
    int k = s % steps;
    double dup = minDup + (maxDup - minDup) * double(i) / double(steps);
    double loss = minLoss + (maxLoss - minLoss) * double(j) / double(steps);
    double trans = minTrans + (maxTrans - minTrans) * double(k) / double(steps);
    jointTree.setRates(dup, loss, trans);
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
  int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  ParallelContext::broadcastDouble(bestRank, bestTrans);
  jointTree.setRates(bestDup, bestLoss, bestTrans);
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
  Logger::info << "best rates " << endl;
  Logger::info << "D " << jointTree.getDupRate() << endl;
  Logger::info << "L " << jointTree.getLossRate() << endl;
  Logger::info << "T " << jointTree.getTransferRate() << endl;
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
  Logger::info << "best rates " << endl;
  Logger::info << "D " << jointTree.getDupRate() << endl;
  Logger::info << "L " << jointTree.getLossRate() << endl;
}


void DTLOptimizer::optimizeDLRatesWindow(JointTree &jointTree) {
  Logger::timed << "Start optimizing DL rates" << endl;
  
  double bestLL = numeric_limits<double>::lowest();
  double newLL = 0;
  double bestDup = 0.0;
  double bestLoss = 0.0;
  double minDup = 0.0;
  double maxDup = 10.0;
  double minLoss = 0.0;
  double maxLoss = 10.0;
  int steps = 10;
  double epsilon = 0.001;
  do {
    bestLL = newLL;
    findBestRatesDL(jointTree, minDup, maxDup, minLoss, maxLoss, steps, bestDup, bestLoss, newLL);
    double offsetDup = (maxDup - minDup) / steps;
    double offsetLoss =(maxLoss - minLoss) / steps;
    minDup = max(0.0, bestDup - offsetDup);
    maxDup = bestDup + offsetDup;
    minLoss = max(0.0, bestLoss - offsetLoss);
    maxLoss = bestLoss + offsetLoss;
  } while (fabs(newLL - bestLL) > epsilon);
  Logger::info << " best rates: " << bestDup << " " << bestLoss <<  " " << newLL << endl;
  if  (!isValidLikelihood(newLL)) {
    Logger::error << "Invalid likelihood " << newLL << endl;
    ParallelContext::abort(10);
  }
}
    
void DTLOptimizer::optimizeDTLRatesWindow(JointTree &jointTree) {
  Logger::timed << "Start optimizing DTL rates" << endl;
  double bestLL = numeric_limits<double>::lowest();
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
  int steps = 5;
  double epsilon = 0.01;
  do {
    bestLL = newLL;
    findBestRatesDTL(jointTree, minDup, maxDup, minLoss, maxLoss, minTrans, maxTrans, steps, bestDup, bestLoss, bestTrans, newLL);
    double offsetDup = (maxDup - minDup) / steps;
    double offsetLoss = (maxLoss - minLoss) / steps;
    double offsetTrans = (maxTrans - minTrans) / steps;
    minDup = max(0.0, bestDup - offsetDup);
    maxDup = bestDup + offsetDup;
    minLoss = max(0.0, bestLoss - offsetLoss);
    maxLoss = bestLoss + offsetLoss;
    minTrans = max(0.0, bestTrans - offsetTrans);
    maxTrans = bestTrans + offsetTrans;
  } while (fabs(newLL - bestLL) > epsilon);
  if  (!isValidLikelihood(newLL)) {
    Logger::error << "Invalid likelihood " << newLL << endl;
    ParallelContext::abort(10);
  }
}

struct DTLRates {
  double rates[3];
  double ll;
  
  DTLRates(double d = 0.0, double l = 0.0, double t = 0.0): ll(0.0) {
    rates[0] = d;
    rates[1] = l;
    rates[2] = t;
  }

  void computeLL(JointTree &jointTree) {
    ensureValidity();
    jointTree.setRates(rates[0], rates[1], rates[2]);
    ll = jointTree.computeReconciliationLoglk();
    if (!isValidLikelihood(ll)) {
      ll = -100000000000;
    }
  }
  
  void computeLL(PerCoreGeneTrees &trees, pll_rtree_t *speciesTree, RecModel model) {
    ensureValidity();
    ll = 0;
    for (auto tree: trees.getTrees()) {
      ReconciliationEvaluation evaluation(speciesTree, tree.mapping, model, true);
      evaluation.setRates(rates[0], rates[1], rates[2]);
      ll += evaluation.evaluate(tree.tree);
    }
    ParallelContext::sumDouble(ll);
    if (!isValidLikelihood(ll)) {
      ll = -100000000000;
    }
    Logger::info << ll << endl;
  }
  
  inline void ensureValidity() {
    rates[0] = max(0.0, rates[0]);
    rates[1] = max(0.0, rates[1]);
    rates[2] = max(0.0, rates[2]);
  }

  inline bool operator <(const DTLRates& v) const {
    return ll > v.ll;
  }
  
  inline bool operator <=(const DTLRates& v) const {
    return ll >=  v.ll;
  }
  
  inline DTLRates operator+(const DTLRates& v) const {
    return DTLRates(rates[0] + v.rates[0], rates[1] + v.rates[1], rates[2] + v.rates[2]);
  }
  
  inline DTLRates operator-(const DTLRates& v) const {
    return DTLRates(rates[0] - v.rates[0], rates[1] - v.rates[1], rates[2] - v.rates[2]);
  }
  
  inline DTLRates operator*(double v) const {
    return DTLRates(v * rates[0], v * rates[1], v * rates[2]);
  }
  
  inline DTLRates operator/(double v) const {
    return DTLRates(rates[0] / v, rates[1] / v, rates[2] / v);
  }
  
  friend ostream& operator<<(ostream& os, const DTLRates &v) {
    os << "(" << v.rates[0] << ", " << v.rates[1] << ", " << v.rates[2] << ", " <<  v.ll  << ")";
    return os;
  }

  inline double distance(const DTLRates v) const {
    double d = 0.0;
    d += pow(rates[0] - v.rates[0], 2.0);
    d += pow(rates[1] - v.rates[1], 2.0);
    d += pow(rates[2] - v.rates[2], 2.0);
    return sqrt(d);
  }
};


/*
 *  Find the point between  r1 and r2. Parallelized over the iterations
 */
DTLRates findBestPoint(DTLRates r1, DTLRates r2, int iterations, JointTree &jointTree) 
{
  DTLRates best = r1;
  best.ll = -100000000000;
  int bestI = 0;
  for (int i = ParallelContext::getRank(); i < iterations; i += ParallelContext::getSize()) {
    DTLRates current = r1 + ((r2 - r1) * (double(i) / double(iterations - 1)));
    current.computeLL(jointTree);
    if (current < best) {
      best = current;
      bestI = i;
    }
  }
  int bestRank = 0;
  ParallelContext::getMax(best.ll, bestRank);
  ParallelContext::broadcastInt(bestRank, bestI);
  if (ParallelContext::getRank() != bestRank) {
    best = r1 + ((r2 - r1) * (double(bestI) / double(iterations - 1)));
    best.ensureValidity();
  }
  ParallelContext::broadcastDouble(bestRank, best.ll);
  return best;
}



void DTLOptimizer::optimizeRateSimplex(JointTree &jointTree, bool transfers)
{
  Logger::timed << "Starting DTL rates optimization" << endl;
  vector<DTLRates> rates;
  rates.push_back(DTLRates(0.01, 0.01, 0.01));
  rates.push_back(DTLRates(1.0, 0.01, 0.01));
  rates.push_back(DTLRates(0.01, 1.0, 1.0));
  if (transfers) {
    rates.push_back(DTLRates(0.01, 0.01, 1.0));
  }
  for (auto &r: rates) {
    r.computeLL(jointTree);
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
    DTLRates xr = findBestPoint(x1, x2, iterations, jointTree);
    if (xr < rates[rates.size() - 1] ) {
      rates.back() = xr;
    }
    currentIt++;
  }
  Logger::timed << "Simplex converged after " << currentIt << " iterations" << endl;
  sort(rates.begin(), rates.end());
  rates[0].computeLL(jointTree);
}
  
void DTLOptimizer::optimizeDTLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model)
{
  Logger::timed << "Todo: merge with other implementations" << endl;
  vector<DTLRates> rates;
  rates.push_back(DTLRates(0.01, 0.01, 0.0));
  rates.push_back(DTLRates(1.0, 0.01, 0.0));
  if (Enums::accountsForTransfers(model)) {
    rates.push_back(DTLRates(0.01, 1.0, 1.0));
    rates.push_back(DTLRates(0.01, 0.01, 1.0));
  } else {
      rates.push_back(DTLRates(0.01, 1.0, 0.0));
  }
  for (auto &r: rates) {
    r.computeLL(geneTrees, speciesTree, model);
  }
  DTLRates worstRate;
  int currentIt = 0;
  while (worstRate.distance(rates.back()) > 0.005) {
    Logger::info << "begin of loop" << endl;
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
    for (int i = 0; i < iterations; i++) {
      DTLRates current = x1 + ((x2 - x1) * (double(i) / double(iterations - 1)));
      current.computeLL(geneTrees, speciesTree, model);
      if (current < bestRates) {
        bestRates = current;
      }
    }
    if (bestRates < rates[rates.size() - 1] ) {
      rates.back() = bestRates;
    }
    currentIt++;
    Logger::info << "end of loop" << endl;
  }
  Logger::timed << "Simplex converged after " << currentIt << " iterations: " << rates[0] << endl;
  sort(rates.begin(), rates.end());
  rates[0].computeLL(geneTrees, speciesTree, model);

}

