#pragma once
#include <random>

class Random {
public:
  Random() = delete;
  static void setSeed(unsigned int seed);
  static int getInt(); 
  static double getProba();

private:
  static std::mt19937_64 _rng;
  static std::uniform_int_distribution<int> _unii;
  static std::uniform_real_distribution<double> _uniproba;
};
