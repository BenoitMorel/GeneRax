#pragma once


class AverageStream {
public:
  AverageStream(unsigned int significantCount = 20):
    _count(0),
    _significantCount(significantCount),
    _average(0)
  {}
  void addValue(double value) {
    _count++;
    _average += (value - _average) / _count;
  }
  double getAverage() const {
    return _average;
  }
  bool isSignificant() const {
    return _count > _significantCount;
  }
private:
  unsigned int _count;
  const unsigned int _significantCount;
  double _average;


};
