#pragma once

#include <string>
#include <IO/Families.hpp>
#include <util/types.hpp>

#include <vector>
#include <unordered_set>

class PLLUnrootedTree;
class UNNIMove;

class MiniBME {
public:
  MiniBME(const PLLUnrootedTree &speciesTree, 
      const Families &families,
      bool pruneMode = true);
  double computeBME(const PLLUnrootedTree &speciesTree);
  double computeNNIDiff(const PLLUnrootedTree &speciesTree,
      const UNNIMove &nni);
private:
  std::vector<DistanceMatrix> _geneDistanceMatrices;
  std::vector<DistanceMatrix> _geneDistanceDenominators;
  std::vector<std::string> _speciesIdToSpeciesString;
  Families _perCoreFamilies;
  std::vector<std::unordered_set<std::string> > _perFamilyCoverageStr;
  std::vector<std::vector<bool> > _perFamilyCoverage;
  StringToUint _speciesStringToSpeciesId;
  bool _prune;
  std::vector<DistanceMatrix> _prunedSpeciesMatrices;
  std::vector<DistanceMatrix> _subBMEs;
  double _computeBMEPrune(const PLLUnrootedTree &speciesTree);
  
  // could be made faster by skipping intersecting subtrees
  // O(n^2)
  void _computeSubBMEs(const PLLUnrootedTree &speciesTree);
  void _computeSubBMEsPrune(const PLLUnrootedTree &speciesTree);
   
};
