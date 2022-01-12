#pragma once

#include <string>
#include <IO/Families.hpp>
#include <util/types.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <vector>
#include <unordered_set>

#include <IO/Logger.hpp>
#include <DistanceMethods/Astrid.hpp>


class Asteroid: public BMEEvaluator {
public:
  Asteroid(const PLLUnrootedTree &speciesTree, 
      const Families &families,
      double minbl);
  virtual ~Asteroid() {}
  /*
   *  Computes the BME score of the species tree
   *  and update the _subBMEs used to speedup the
   *  getBestSPR call.
   */
  virtual double computeBME(const PLLUnrootedTree &speciesTree);

  /**
   *  Tries all possible SPR moves from the speciesTree
   *  and returns the best one as well as the difference
   *  between the current score and the new score
   */
  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves);

private:
  // _getBestSPRRecMissing[k][i][j] is the distance
  // between species i and j for the gene family k
  // this value is only valid if 
  // _geneDistanceDenominators[k][i][j] != 0.0
  std::vector<DistanceMatrix> _geneDistanceMatrices;
  // see _geneDistanceMatrices
  std::vector<DistanceMatrix> _geneDistanceDenominators;
  // number of species nodes
  size_t _N;
  // number of families 
  size_t _K;
  // _perFamilyCoverageStr[k] is the set of labels of the species
  // covered by the family k
  std::vector<std::unordered_set<std::string> > _perFamilyCoverageStr;
  // _perFamilyCoverage[k][i] is true if the family k covers the species i
  std::vector<std::vector<bool> > _perFamilyCoverage;
  // species label to species ID (which is NOT the node_index)
  StringToUint _speciesStringToSpeciesId;
  // _prunedSpeciesMatrices[k] is the internode distance for the
  // species tree induced by the family k
  std::vector<DistanceMatrix> _prunedSpeciesMatrices;
  // _subBMEs[k][i][j] == average distance between the subtrees
  // i and j (indexed with node_index). Only defined for non-
  // intersecting i and j.

#ifdef OLDBME
  std::vector<DistanceMatrix> _subBMEs;
  double getCell(unsigned int sp1, unsigned int sp2, unsigned int k) {
    return _subBMEs[sp1][sp2][k];
  }
  void setCell(unsigned int sp1, unsigned int sp2, unsigned int k, double v) {
    _subBMEs[sp1][sp2][k] = v;
  }
#else
  std::vector<double> _subBMEs;
  double getCell(size_t sp1, size_t sp2, size_t k) {
    size_t index = ((sp1 * _N) + sp2) * _K + k;
    return _subBMEs[index];
  }
  void setCell(size_t sp1, size_t sp2, size_t k, double v) {
    size_t index = ((sp1 * _N) + sp2) * _K + k;
    _subBMEs[index] = v;
  }
#endif
  
  // _hasChildren[i][k] == true if there is at least one leaf
  // under i that belongs to the species tree induced by the
  // family k
  BoolMatrix _hasChildren;
  // _belongsToPruned[i][k] == true if the node i belongs to
  // the species tree induced by the family k
  BoolMatrix _belongsToPruned;
  // _pows[i] == pow(2, i) (precomputed to speedup computations)
  std::vector<double> _pows;

  struct SubBMEToUpdate {
    corax_unode_t *pruned;
    corax_unode_t *before;
    corax_unode_t *after;
    std::vector<corax_unode_t *> between;
    std::vector<corax_unode_t *> tempBetween;
    SubBMEToUpdate() {reset();}
    void reset() {
      pruned = nullptr;
      before = nullptr;
      after = nullptr;
      tempBetween.clear();
      between.clear();
    }
  };
  SubBMEToUpdate _toUpdate;

  double _computeBMEPrune(const PLLUnrootedTree &speciesTree);
  
  // O(n^2)
  void _computeSubBMEsPrune(const PLLUnrootedTree &speciesTree);

  struct StopCriterion {
    unsigned int maxRadius;
    unsigned int maxRadiusWithoutImprovement;
    unsigned int noImprovement;
    StopCriterion(): maxRadius(99999),
      maxRadiusWithoutImprovement(4),
      noImprovement(0)
      {}
  };

  // lot of copies hgre...
  void _getBestSPRRecMissing(unsigned int s,
     StopCriterion stopCriterion,
     std::vector<unsigned int> sprime, 
     std::vector<corax_unode_t *> W0s, 
     corax_unode_t *Wp, 
     corax_unode_t *Wsminus1, 
     corax_unode_t *Vsminus1, 
     std::vector<double> delta_Vsminus2_Wp, // previous deltaAB
     corax_unode_t *Vs, 
     double Lsminus1, // L_s-1
     corax_unode_t *&bestRegraftNode,
     SubBMEToUpdate &subBMEToUpdate,
     double &bestLs,
     unsigned int &bestS,
     const BoolMatrix &belongsToPruned,
     const BoolMatrix &hasChildren,
     std::vector<bool> Vsminus2HasChildren, // does Vsminus2 have children after the previous moves
     std::vector<bool> Vsminus1HasChildren); // does Vsminus1 have children after the previous moves
   

  // computes _subBMEs[k][i1][i2] for all k
  void _computeSubBMEsPruneRec(corax_unode_t *n1,
    corax_unode_t *n2,
    BoolMatrix &treated);
  bool  getBestSPRFromPrune(unsigned int maxRadiusWithoutImprovement,
      corax_unode_t *prunedNode,
      corax_unode_t *&bestRegraftNode,
      double &bestDiff,
      unsigned int &bestS);
};
