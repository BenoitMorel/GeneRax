#pragma once

#include <string>
#include <IO/Families.hpp>

class PLLUnrootedTree;

class BME {
  public:
  static double computeBME(const PLLUnrootedTree &speciesTree,
      const Families &families);

};
