#pragma once

#include <vector>
#include <string>
#include <IO/ParallelOfstream.hpp>
#include <IO/Logger.hpp>
#include <fstream>
#include <cmath>


class Parameters {
public:
  Parameters(): _score(0.0) {}
  
  Parameters(unsigned int dimensions): _parameters(dimensions, 0.0),
    _score(0.0)
  {
  }

  Parameters(unsigned int number, const Parameters &initValue):
    _score(0.0)
  {
    _parameters.resize(number * initValue.dimensions());
    unsigned int index = 0;
    for (unsigned int i = 0; i < number; ++i) {
      for (auto p: initValue._parameters) {
        _parameters[index++] = p;
      }
    }
  }
  
  Parameters(double d, double l): _score(0.0) {
    _parameters.push_back(d);
    _parameters.push_back(l);
  }
  
  Parameters(double d, double l, double t): _score(0.0) {
    _parameters.push_back(d);
    _parameters.push_back(l);
    _parameters.push_back(t);
  }
  
  Parameters(double d, double l, double t, double i): _score(0.0) {
    _parameters.push_back(d);
    _parameters.push_back(l);
    _parameters.push_back(t);
    _parameters.push_back(i);
  }

  Parameters(const std::vector<double> &parameters): _parameters(parameters),
    _score(0.0)
  {
  }

  Parameters(const std::string &src): _score(0.0) {load(src);} 
  
  inline unsigned int dimensions() const {
    return static_cast<unsigned int>(_parameters.size());
  }

  inline void addValue(double value) {
    _parameters.push_back(value);
  }

  inline Parameters getSubParameters(unsigned int begin, unsigned int size) const {
    return Parameters(std::vector<double>(_parameters.begin() + begin, 
          _parameters.begin() + begin + size));
  }

  void ensurePositivity() {
    for (auto &p: _parameters) {
      p = std::max(0.0000001, p);
      //p = std::min(1.0, p);
    }
  }

  double operator [](unsigned int i) const    {return _parameters[i];}
  double & operator [](unsigned int i) {return _parameters[i];}
  
  inline double getScore() const {return _score;}
  inline void setScore(double score) {_score = score;}
  
  inline bool operator <(const Parameters& v) const {
    return getScore() > v.getScore();
  }
  
  inline bool operator <=(const Parameters& v) const {
    return getScore() >=  v.getScore();
  }
  
  inline Parameters operator+(const Parameters& v) const {
    assert(dimensions() == v.dimensions());
    Parameters res(*this);
    for (unsigned int i = 0; i < dimensions(); ++i) {
      res._parameters[i] += v._parameters[i];
    }
    return res;
  }
  
  inline Parameters operator-(const Parameters& v) const {
    assert(dimensions() == v.dimensions());
    Parameters res(*this);
    for (unsigned int i = 0; i < dimensions(); ++i) {
      res._parameters[i] -= v._parameters[i];
    }
    return res;
  }
  
  inline Parameters operator*(double v) const {
    Parameters res(*this);
    for (unsigned int i = 0; i < dimensions(); ++i) {
      res._parameters[i] *= v;
    }
    return res;
  }
  
  inline Parameters operator/(double v) const {
    Parameters res(*this);
    for (unsigned int i = 0; i < dimensions(); ++i) {
      res._parameters[i] /= v;
    }
    return res;
  }
  
  inline double distance(const Parameters &v) const {
    double d = 0.0;
    for (unsigned int i = 0; i < dimensions(); ++i) {
      d += std::pow((*this)[i] - v[i], 2.0);
    }
    return std::sqrt(d);
  }

  void normalize(double norm = 1.0) {
    double av = distance(Parameters(dimensions()));
    *this = *this * (norm / av);
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Parameters &v) {
    os << "(";
    for (auto value: v._parameters) {
      os << value << ", ";
    }
    os << "score = " << v.getScore() << ")";
    return os;
  }
  
  void save(const std::string &dest) 
  {
    ParallelOfstream os(dest);
    for (auto &r: _parameters) {
      os << r << " ";
    }
  }

  void load(const std::string &src) 
  {
    _parameters.clear();
    std::ifstream is2(src);
    std::string debugStr;
    std::getline(is2, debugStr);
    std::ifstream is(src);
    if (!is.good()) {
      std::cerr << "Cannot read parameters from " << src << std::endl;
    }
    assert(is.good());
    while (!is.eof()) {
      double a = -42;
      is >> a;
      if (!is.good())  {
        break;
      }
      _parameters.push_back(a);
    }
  }

  const std::vector<double> &getVector() const {return _parameters;}
private:
  std::vector<double> _parameters;
  double _score;
};



