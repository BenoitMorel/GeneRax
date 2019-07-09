#pragma once

#include <cassert>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <IO/ParallelOfstream.hpp>
#include <IO/Logger.hpp>
#include <stdio.h>
#include <stdlib.h>

struct DTLRates {
  double rates[3];
  double ll;
  
  DTLRates(double d = 0.0, double l = 0.0, double t = 0.0): ll(0.0) {
    rates[0] = d;
    rates[1] = l;
    rates[2] = t;
  }


  inline void ensureValidity() {
    rates[0] = std::max(0.0, rates[0]);
    rates[1] = std::max(0.0, rates[1]);
    rates[2] = std::max(0.0, rates[2]);
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
  
  friend std::ostream& operator<<(std::ostream& os, const DTLRates &v) {
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

static double randfrom(double min, double max, 
    std::default_random_engine &gen)
{
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  if (distribution(gen) < 0.5) {
    return min;
  } else {
    return max;
  }
}

class DTLRatesVector {
public:
  DTLRatesVector(): _ll(0.0) {}
  DTLRatesVector(unsigned int size): _rates(size), _ll(0.0) {}
  DTLRatesVector(unsigned int size, const DTLRates &rate): _rates(size, rate), _ll(0.0) {}  
  DTLRatesVector(const std::string &src): _ll(0.0) {load(src);} 
  DTLRatesVector(const DTLRates &rates): _ll(0.0) {_rates.push_back(rates);}
  void initRandom(double min, double max, std::default_random_engine &gen) {
    for (auto &rate: _rates) {
      rate = DTLRates(randfrom(min, max, gen), randfrom(min, max, gen), randfrom(min, max, gen));
    }
  }
  
  inline bool operator <(const DTLRatesVector& v) const {
    return _ll > v._ll;
  }
  
  inline bool operator <=(const DTLRatesVector& v) const {
    return _ll >=  v._ll;
  }
 
  inline double getLL() const {return _ll;}
  
  inline void setLL(double ll) { _ll = ll;}

  inline unsigned int size() const {return _rates.size();}

  inline const DTLRates &getRates(unsigned int index) const {return _rates[index];}
  
  inline DTLRates &getRates(unsigned int index) {return _rates[index];}

  inline const std::vector<DTLRates> &getRatesVector() const {return _rates;}

  inline void ensureValidity() {
    for (auto &r: _rates) {
      r.ensureValidity();
    }
  }

  inline DTLRatesVector operator+(const DTLRatesVector& v) const {
    assert(size() == v.size());
    DTLRatesVector res(size());
    for (unsigned int i = 0; i < size(); ++i) {
      res.getRates(i) = getRates(i) + v.getRates(i); 
    }
    return res;
  }
  
  inline DTLRatesVector operator-(const DTLRatesVector& v) const {
    assert(size() == v.size());
    DTLRatesVector res(size());
    for (unsigned int i = 0; i < size(); ++i) {
      res.getRates(i) = getRates(i) - v.getRates(i); 
    }
    return res;
  }
  
  inline DTLRatesVector operator*(double v) const {
    DTLRatesVector res(size());
    for (unsigned int i = 0; i < size(); ++i) {
      res.getRates(i) = getRates(i) * v; 
    }
    return res;
  }
  
  inline DTLRatesVector operator/(double v) const {
    DTLRatesVector res(size());
    for (unsigned int i = 0; i < size(); ++i) {
      res.getRates(i) = getRates(i) / v; 
    }
    return res;
  }
 
  friend std::ostream& operator<<(std::ostream& os, const DTLRatesVector &v) {
    os << "[";
    for (auto &rates: v._rates) {
      os << rates << " ";
    }
    return os;
  }

  void save(const std::string &dest) 
  {
    ParallelOfstream os(dest);
    for (auto &r: _rates) {
      os << r.rates[0] << " " << r.rates[1] << " " << r.rates[2] << std::endl;
    }
  }

  void load(const std::string &src) 
  {
    _rates.clear();
    std::ifstream is(src);
    std::string line;
    while (std::getline(is, line))
    {
      if (line.size() == 0) {
        break;
      }
      std::istringstream iss(line);
      double a, b, c;
      iss >> a >> b >> c;
      _rates.push_back(DTLRates(a, b, c));
    }
  }
  
  inline double distance(const DTLRatesVector &v) const {
    assert(v.size() == size());
    double res = 0.0;
    for (unsigned int i = 0; i < size(); ++i) {
      res += pow(getRates(i).distance(v.getRates(i)), 2.0);
    }
    return sqrt(res) / double(size());
  }

private:
  std::vector<DTLRates> _rates;
  double _ll;
};
