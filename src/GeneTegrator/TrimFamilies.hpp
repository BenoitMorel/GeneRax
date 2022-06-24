#pragma once

#include <IO/Families.hpp>

class TrimFamilies {
public: 
  static void trimHighCladesNumber(Families &families,
      double keepRatio);

  static void trimMinSpeciesCoverage(Families &families,
      unsigned int minCoverage);
};


