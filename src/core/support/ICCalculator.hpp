#pragma once
#include <vector>
#include <string>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <util/types.hpp>

class ICCalculator {
public:
  ICCalculator(const std::string &referenceTreePath,
      const Families &families,
      bool paralogy);

};
