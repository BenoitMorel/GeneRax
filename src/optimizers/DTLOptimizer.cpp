#include <optimizers/DTLOptimizer.hpp>

#include <ParallelContext.hpp>
#include <treeSearch/JointTree.hpp>
#include <IO/Logger.hpp>

#include <limits>
#include <algorithm>

bool isValidLikelihood(double ll) {
  return isnormal(ll) && ll < -0.0000001;
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



DTLRates findBestPoint(DTLRates r1, DTLRates r2, JointTree &jointTree) 
{
  //r1.computeLL(jointTree);
  DTLRates best = r1;
  best.ll = -100000000000;
  int iterations = 10;
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


void DTLOptimizer::optimizeDTLRates(JointTree &jointTree, bool transfers)
{
  Logger::timed << "Starting DTL rates optimization" << endl;
  vector<DTLRates> rates;
  rates.push_back(DTLRates(0.01, 0.01, 0.0));
  rates.push_back(DTLRates(1.0, 0.01, 0.0));
  rates.push_back(DTLRates(0.01, 1.0, 0.0));
  if (transfers) {
    rates.push_back(DTLRates(0.01, 0.01, 1.0));
  }
  for (auto &r: rates) {
    r.computeLL(jointTree);
  }
  double alpha = 1.0;
  double gamma = 2.0;
  double rho = 0.5;
  double sigma = 0.5;


  while (rates[0].distance(rates.back()) > 0.001) {
    sort(rates.begin(), rates.end());
    for (auto &rate: rates)  {
      Logger::info << rate << " ";
    }
    Logger::info << endl;
    // centroid
    DTLRates x0;
    for (unsigned int i = 0; i < rates.size() - 1; ++i) {
      x0 = x0 + rates[i];
    }
    x0 = x0 / double(rates.size() - 1);
    // reflexion
    DTLRates xr = x0 + (x0 - rates.back()) * alpha;  
    xr = findBestPoint(x0, xr, jointTree);
    if (rates[0] <= xr && xr < rates[rates.size() - 2] ) {
      Logger::info << "reflexion" << endl;
      rates.back() = xr;
      continue;
    }
    // expansion
    if (xr < rates[0]) {
      DTLRates xe = x0 + (xr - x0) * gamma;
      xe = findBestPoint(xe, xr, jointTree);
      //xe = findBestRatesDL(x0, xr, jointTree);
      //xe.computeLL(jointTree);
      if (xe < xr) {
        Logger::info << "expansion" << endl;
        rates.back() = xe;
      } else {
        Logger::info << "reflexion/expansion" << endl;
        rates.back() = xr;
      }
      continue;
    }
    // contraction
    DTLRates xc = x0 + (rates.back() - x0) * rho;
    xc.computeLL(jointTree); 
    if (xc < rates.back()) {
      rates.back() = xc;
      Logger::info << "contraction" << endl;
      continue;
    }
    // shrink
    Logger::info << "shrink" << endl;
    for (unsigned int r = 1; r < rates.size(); ++r) {
      rates[r] = rates[0] + (rates[r] - rates[0]) * sigma;
      rates[r].computeLL(jointTree);
    }
  }
  sort(rates.begin(), rates.end());
  rates[0].computeLL(jointTree);
  Logger::info << " best rates: " << rates[0].rates[0] << " " << rates[0].rates[1] << " " << rates[0].rates[2] << " " << rates[0].ll <<  endl;
}

