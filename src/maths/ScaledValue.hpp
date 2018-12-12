#ifndef __SCALED_VALUE_HPP__ 
#define __SCALED_VALUE_HPP__

#include <limits>
#include <climits>
#define JS_SCALE_FACTOR 115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define JS_SCALE_THRESHOLD (1.0/JS_SCALE_FACTOR)
#include <iostream>


/*
 *  Scaled double, for operations on very small values
 */
class ScaledValue {
public:
  
  ScaledValue():value(0), scaler(INT_MAX) {
  }

  ScaledValue(double v, int s):value(v), scaler(s) {
    if (value == 0) {
      scaler = INT_MAX;
      value = 0;
    }
    else {
      while (value < JS_SCALE_THRESHOLD) {
        scaler += 1;
        value *= JS_SCALE_FACTOR;
      }
    }
  } 
 
  ScaledValue& operator/=(double v) {
    value /= v;
  }

  ScaledValue& operator+=(const ScaledValue& v) {
    if (v.scaler == scaler) {
      value += v.value;
    } else if (v.scaler < scaler) {
      value = v.value;
      scaler = v.scaler;
    }
    return *this;
  }
  
  ScaledValue operator+(const ScaledValue& v) {
    if (v.scaler == scaler) {
      return ScaledValue(v.value + value, scaler);
    } else if (v.scaler < scaler) {
      return v;
    } else {
      return *this;
    }
  }
  
  ScaledValue operator*(const ScaledValue& v) {
    return ScaledValue (v.value * value, v.scaler + scaler);  
  }

  ScaledValue operator*(double v) {
    return ScaledValue(v * value, scaler);
  }

  bool operator <(const ScaledValue& v)
  {
    if (scaler != v.scaler) {
      return scaler > v.scaler;
    }
    return value < v.value;
  }

  double getLogValue() {
    if (scaler == INT_MAX) {
      return -std::numeric_limits<double>::infinity();
    }
    return log(value) + scaler * log(JS_SCALE_THRESHOLD);
  }

  friend ostream& operator<<(ostream& os, const ScaledValue &v) {
    os << "(" << v.value << "," << v.scaler << ")";
    return os;
  }

private:
  double value;
  int scaler;
};



#endif
