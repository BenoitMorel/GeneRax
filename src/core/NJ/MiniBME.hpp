#pragma once

#include <string>
#include <IO/Families.hpp>
#include <util/types.hpp>

#include <vector>
#include <unordered_set>

class PLLUnrootedTree;

class MiniBME {
public:
  MiniBME(const PLLUnrootedTree &speciesTree, 
      const Families &families,
      bool pruneMode = true);
  double computeBME(const PLLUnrootedTree &speciesTree);

private:
  std::vector<DistanceMatrix> _geneDistanceMatrices;
  std::vector<DistanceMatrix> _geneDistanceDenominators;
  std::vector<std::string> _speciesIdToSpeciesString;
  Families _perCoreFamilies;
  std::vector<std::unordered_set<std::string> > _perFamilyCoverage;
  StringToUint _speciesStringToSpeciesId;
  bool _prune;
  double _computeBMEPrune(const PLLUnrootedTree &speciesTree);
};
