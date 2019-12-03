#pragma once

#include <limits>
#include <climits>
#include <iostream>

#define JS_SCALE_FACTOR \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define JS_SCALE_THRESHOLD (1.0/JS_SCALE_FACTOR)

const int NULL_SCALER = INT_MAX / 2 - 1;

/*
 *  Class representing a double value with a high precision. It stores a double
 *  value, and a integer scaling value to represent very small values
 *
 *  When the value is null, the scaler is set to NULL_SCALER
 */
class ScaledValue {
public:
  /**
   * Constructor for a null value
   */
  ScaledValue():value(0.0), scaler(NULL_SCALER) {
  }

  void setNull() {
    value = 0.0;
    scaler = NULL_SCALER;
  }

  void checkNull() {
    if (value == 0.0) {
      scaler = NULL_SCALER;
    }
  }

  void scale() {
    if (value < JS_SCALE_THRESHOLD) {
      scaler += 1;
      value *= JS_SCALE_FACTOR;
    }
  }
  
  explicit ScaledValue(double v):value(v), scaler(0)  {
    checkNull();
  } 

  /**
   * General constructor
   * @param v value
   * @param s scaler
   */
  explicit ScaledValue(double v, int s):value(v), scaler(s)  {
    checkNull();
  } 


  /**
   * ScaledValue sum operator
   */
  inline ScaledValue& operator+=(const ScaledValue& v) {
    if (v.scaler == scaler) {
      value += v.value;
    } else if (v.scaler < scaler) {
      value = v.value;
      scaler = v.scaler;
    }
    scale();
    return *this;
  }

  
  /**
   * ScaledValue sum operator
   */
  inline ScaledValue operator+(const ScaledValue& v) const {
    if (v.scaler == scaler) {
      return ScaledValue(v.value + value, scaler);
    } else if (v.scaler < scaler) {
      return v;
    } else {
      return *this;
    }
  }
  
  /**
   * ScaledValue minux operator
   */
  inline ScaledValue operator-(const ScaledValue& v) const {
    if (v.scaler == scaler) {
      if (value - v.value < 0.0) {
        if (fabs(value - v.value) < 0.0000000001) {
          return ScaledValue();
        }
        std::cerr.precision(17);
        std::cerr << *this << " - " << v << std::endl;
      }
      assert(value - v.value >= 0);
      auto res = ScaledValue(value - v.value, scaler);
      res.scale();
      return res;
    } else if (v.scaler < scaler) {
      std::cerr << *this << " - " << v << std::endl;
      assert(false); // we do not allow negative values for now
      return v;
    } else {
      return *this;
    }
  }
  
  /**
   * ScaledValue multiplication operator
   */
  inline ScaledValue operator*(const ScaledValue& v) const {
    auto res = ScaledValue (v.value * value, v.scaler + scaler);
    res.scale();
    res.checkNull();
    return res;
  }
  
  /**
   * ScaledValue multiplication operator
   */
  inline ScaledValue& operator*=(const ScaledValue& v) {
    value *= v.value;
    if (scaler != NULL_SCALER && v.scaler != NULL_SCALER) {
      scaler += v.scaler;   
      scale();
    } else {
      setNull();
    }
    return *this;
  }

  /**
   * double multiplication operator
   */
  inline ScaledValue operator*(double v) const {
    auto res = ScaledValue(v * value, scaler);
    res.checkNull();
    res.scale();
    return res;
  }

  /**
   * double multiplication operator
   */
  inline ScaledValue& operator*=(double v) {
    value *= v;
    checkNull();
    scale();
    return *this;
  }
  
  /**
   * Division operator
   */
  inline ScaledValue& operator/=(double v) {
    value /= v; 
    return *this;
  }


  /**
   * Comparison with ScaledValue operator
   */
  inline bool operator <(const ScaledValue& v) const
  {
    if (isNull()) {
      return !v.isNull();
    }
    if (scaler != v.scaler) {
      return scaler > v.scaler;
    }
    return value < v.value;
  }
  
  inline bool operator ==(const ScaledValue& v) const
  {
    return scaler == v.scaler && (std::fabs(v.value-value) <= std::numeric_limits<double>::epsilon());
  }

  inline bool operator <=(const ScaledValue& v) const
  {
    if (isNull()) {
      return true;
    }
    if (scaler != v.scaler) {
      return scaler > v.scaler;
    }
    return value <= v.value;
  }

  /**
   * @return true if value is 0
   */
  inline bool isNull() const {
    return value == 0.0;
  }

  /**
   * @return the logarithm value as a double
   */
  inline double getLogValue() {
    if (scaler == NULL_SCALER) {
      return -std::numeric_limits<double>::infinity();
    }
    return log(value) + scaler * log(JS_SCALE_THRESHOLD);
  }

  
  bool isProba() const {
    return *this <= ScaledValue(1.0) && ScaledValue() <= *this;
  }

  /**
   *  std::ofstream operator
   */
  friend std::ostream& operator<<(std::ostream& os, const ScaledValue &v) {
    os << "(" << v.value << "," << v.scaler << ")";
    return os;
  }
 

  double value;
  int scaler;
};


