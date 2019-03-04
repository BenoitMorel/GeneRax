#include <optimizers/DTLOptimizer.hpp>

#include <ParallelContext.hpp>
#include <treeSearch/JointTree.hpp>
#include <IO/Logger.hpp>

#include <limits>
#include <algorithm>

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


void DTLOptimizer::optimizeDTLRates(JointTree &jointTree, const string &method) {
  if (method == "window") {
    optimizeDTLRatesWindow(jointTree);
  } else if (method == "simplex") {
    optimizeRateSimplex(jointTree, true);
  } else {
    assert(false);
  }
}

void DTLOptimizer::optimizeDLRates(JointTree &jointTree, const string &method) {
  if (method == "window") {
    optimizeDLRatesWindow(jointTree);
  } else if (method == "simplex") {
    optimizeRateSimplex(jointTree, false);
  } else {
    assert(false);
  }
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
    rates[0] = max(0.0, rates[0]);
    rates[1] = max(0.0, rates[1]);
    rates[2] = max(0.0, rates[2]);
    jointTree.setRates(rates[0], rates[1], rates[2]);
    ll = jointTree.computeReconciliationLoglk();
    if (!isValidLikelihood(ll)) {
      ll = -100000000000;
    }

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



DTLRates findBestPoint(DTLRates r1, DTLRates r2, JointTree &jointTree) 
{
  r1.computeLL(jointTree);
  DTLRates best = r1;
  int iterations = 10;
  for (int i = 0; i < iterations; ++i) {
    double ratio = double(i) / double(iterations);
    DTLRates current = r1 + ((r2 - r1) * ratio);
    current.computeLL(jointTree);
    //Logger::info << best.ll << " " << current.ll << endl;
    if (current < best) {
      best = current;
    }
  }
  return best;
}


void DTLOptimizer::optimizeRateSimplex(JointTree &jointTree, bool transfers)
{
  vector<DTLRates> rates;
  rates.push_back(DTLRates(0.01, 0.01));
  rates.push_back(DTLRates(1.0, 0.01));
  rates.push_back(DTLRates(0.01, 1.0));
  if (transfers) {
    rates.push_back(DTLRates(0.01, 0.1, 1.0));
  }
  for (auto &r: rates) {
    r.computeLL(jointTree);
  }
  double alpha = 1.0;
  double gamma = 2.0;
  double rho = 0.5;
  double sigma = 0.5;


  while (rates[0].distance(rates.back()) > 0.05) {
  //for (int i = 0; i < 10; ++i) {
    sort(rates.begin(), rates.end());
    // centroid
    DTLRates x0;
    for (int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion
    DTLRates xr = x0 + (x0 - rates.back()) * alpha;  
    xr.computeLL(jointTree);
    if (rates[0] <= xr && xr < rates[rates.size() - 2] ) {
      rates.back() = xr;
      continue;
    }
    // expansion
    if (xr < rates[0]) {
      DTLRates xe = x0 + (xr - x0) * gamma;
      //xe = findBestRatesDL(x0, xr, jointTree);
      xe.computeLL(jointTree);
      if (xe < xr) {
        rates.back() = xe;
      } else {
        rates.back() = xr;
      }
      continue;
    }
    // contraction
    DTLRates xc = x0 + (rates.back() - x0) * rho;
    xc.computeLL(jointTree); 
    if (xc < rates.back()) {
      rates.back() = xc;
      continue;
    }
    // shrink
    for (unsigned int r = 1; r < rates.size(); ++r) {
      rates[r] = rates[0] + (rates[r] - rates[0]) * sigma;
      rates[r].computeLL(jointTree);
    }
  }
  sort(rates.begin(), rates.end());
  rates[0].computeLL(jointTree);
}

