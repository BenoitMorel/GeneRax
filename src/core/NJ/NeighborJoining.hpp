#pragma once

#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>
#include <vector>
#include <string>
#include <util/types.hpp>

class NeighborJoining {
public:
  static std::unique_ptr<PLLRootedTree> applyNJ(
      DistanceMatrix &distanceMatrix,
      std::vector<std::string> &speciesIdToSpeciesString,
      StringToUint &speciesStringToSpeciesId);

};


