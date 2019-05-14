#pragma once

struct DTLRates {
  double rates[3];
  double ll;
  
  DTLRates(double d = 0.0, double l = 0.0, double t = 0.0): ll(0.0) {
    rates[0] = d;
    rates[1] = l;
    rates[2] = t;
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

