#pragma once

#include <random>

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

class DTLRatesVector {
public:
  DTLRatesVector(unsigned int size): _rates(size), _ll(0.0) {}
  DTLRatesVector(unsigned int size, const DTLRates &rate): _rates(size, rate), _ll(0.0) {}  
  
  void initRandom(double min = 0.00001, double max = 1.0) {
    std::random_device rd;  
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    for (auto &rate: _rates) {
      rate = DTLRates(dis(gen), dis(gen), dis(gen));
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
  
  inline double distance(const DTLRatesVector &v) const {
    assert(v.size() == size());
    double res = 0.0;
    for (unsigned int i = 0; i < size(); ++i) {
      res += getRates(i).distance(v.getRates(i));
    }
    return res;
  }

private:
  std::vector<DTLRates> _rates;
  double _ll;
};
