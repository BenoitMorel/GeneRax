#pragma once

#include <string>
#include <IO/Families.hpp>
#include <util/types.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <vector>
#include <unordered_set>

using BoolMatrix = std::vector< std::vector<bool> >;

class BMEEvaluator {
public:
  virtual ~BMEEvaluator() {}
  virtual double computeBME(const PLLUnrootedTree &speciesTree) = 0;
  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      corax_unode_t *&bestPruneNode,
      corax_unode_t *&bestRegraftNode,
      double &bestDiff) = 0;
};

class MiniBME: public BMEEvaluator {
public:
  MiniBME(const PLLUnrootedTree &speciesTree, 
      const Families &families);
  virtual ~MiniBME() {}
  virtual double computeBME(const PLLUnrootedTree &speciesTree);
  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      corax_unode_t *&bestPruneNode,
      corax_unode_t *&bestRegraftNode,
      double &bestDiff);
private:
  std::vector<DistanceMatrix> _geneDistanceMatrices;
  std::vector<DistanceMatrix> _geneDistanceDenominators;
  std::vector<std::string> _speciesIdToSpeciesString;
  Families _perCoreFamilies;
  std::vector<std::unordered_set<std::string> > _perFamilyCoverageStr;
  std::vector<std::vector<bool> > _perFamilyCoverage;
  StringToUint _speciesStringToSpeciesId;
  std::vector<DistanceMatrix> _subBMEs;
  BoolMatrix _hasChildren;
  BoolMatrix _belongsToPruned;
  std::vector<double> _pows;
  
  bool  getBestSPRFromPrune(corax_unode_t *prunedNode,
      corax_unode_t *&bestRegraftNode,
      double &bestDiff,
      unsigned int &bestS);
  // could be made faster by skipping intersecting subtrees
  // O(n^2)
  void _computeSubBMEs(const PLLUnrootedTree &speciesTree);
   
};
