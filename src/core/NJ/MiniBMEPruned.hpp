#pragma once

#include <string>
#include <IO/Families.hpp>
#include <util/types.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <vector>
#include <unordered_set>

#include <NJ/MiniBME.hpp>

class MiniBMEPruned: public BMEEvaluator {
public:
  MiniBMEPruned(const PLLUnrootedTree &speciesTree, 
      const Families &families);
  virtual ~MiniBMEPruned() {}
  virtual double computeBME(const PLLUnrootedTree &speciesTree);

  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      pll_unode_t *&bestPruneNode,
      pll_unode_t *&bestRegraftNode,
      double &bestDiff);
private:
  std::vector<DistanceMatrix> _geneDistanceMatrices;
  std::vector<DistanceMatrix> _geneDistanceDenominators;
  std::vector<std::string> _speciesIdToSpeciesString;
  Families _perCoreFamilies;
  std::vector<std::unordered_set<std::string> > _perFamilyCoverageStr;
  std::vector<std::vector<bool> > _perFamilyCoverage;
  StringToUint _speciesStringToSpeciesId;
  std::vector<DistanceMatrix> _prunedSpeciesMatrices;
  std::vector<DistanceMatrix> _subBMEs;
  BoolMatrix _hasChildren;
  BoolMatrix _belongsToPruned;
  std::vector<double> _pows;
  double _computeBMEPrune(const PLLUnrootedTree &speciesTree);
  
  // could be made faster by skipping intersecting subtrees
  // O(n^2)
  void _computeSubBMEsPrune(const PLLUnrootedTree &speciesTree);
 
  void _getBestSPRRecMissing(unsigned int s,
     std::vector<unsigned int> sprime, // copy!!
     std::vector<pll_unode_t *> W0s, 
     pll_unode_t *Wp, 
     pll_unode_t *Wsminus1, 
     pll_unode_t *Vsminus1, 
     std::vector<double> delta_Vsminus2_Wp, // previous deltaAB
     pll_unode_t *Vs, 
     double Lsminus1, // L_s-1
     pll_unode_t *&bestRegraftNode,
     double &bestLs,
     unsigned int &bestS,
     const std::vector<DistanceMatrix> &subBMEs,
     const BoolMatrix &belongsToPruned,
     const BoolMatrix &hasChildren,
     std::vector<bool> Vsminus2HasChildren, // does Vsminus2 have children after the previous moves
     std::vector<bool> Vsminus1HasChildren); // does Vsminus1 have children after the previous moves
   
  bool  getBestSPRFromPrune(pll_unode_t *prunedNode,
      pll_unode_t *&bestRegraftNode,
      double &bestDiff,
      unsigned int &bestS);
};
