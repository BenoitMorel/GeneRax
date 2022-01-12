#pragma once

#include <string>
#include <IO/Families.hpp>
#include <util/types.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <vector>
#include <unordered_set>

using BoolMatrix = std::vector< std::vector<bool> >;

struct SPRMove {
  corax_unode_t *pruneNode;
  corax_unode_t *regraftNode;
  double score;
  SPRMove(corax_unode_t *pruneNode,
      corax_unode_t *regraftNode,
      double score): 
    pruneNode(pruneNode),
    regraftNode(regraftNode),
    score(score) {}
  
  bool operator < (const SPRMove& other) const {
    return other.score < score;
    /*
    if (other.score < score) {
      return true;
    } else if (other.score > score) {
      return false;
    }
    return other.prunedNode->node_index + other.regraftNode->node_index <
      pruneNode->node_index + regraftNode->node_index;
  */
  }
  bool operator ==(const SPRMove& other) const {
    return other.score == score;
  }
};

class BMEEvaluator {
public:
  virtual ~BMEEvaluator() {}
  virtual double computeBME(const PLLUnrootedTree &speciesTree) = 0;
  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves) = 0;
};

class Astrid: public BMEEvaluator {
public:
  Astrid(const PLLUnrootedTree &speciesTree, 
      const Families &families,
      double minbl);
  virtual ~Astrid() {}
  virtual double computeBME(const PLLUnrootedTree &speciesTree);
  virtual void getBestSPR(PLLUnrootedTree &speciesTree,
      unsigned int maxRadiusWithoutImprovement,
      std::vector<SPRMove> &bestMoves);
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
