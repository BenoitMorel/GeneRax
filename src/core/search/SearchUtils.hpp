#pragma once

#include <search/Moves.hpp>

#include <memory>

#include <unordered_map>

class JointTree;

class SearchUtils {
public:
  /**
   *  Apply a given move, compute its likelihood, and rollback the move
   *  Return true if the move improves the likelihood
   *
   *
   *  @param jointTree the current tree
   *  @param move The move to test
   *  @param newLoglk The likelihood of the new tree
   *  @param blo Do we apply branch length optimization?
   *
   *  Parallelization: this function is local to one rank
   */
  static void testMove(JointTree &jointTree,
    SPRMove &move,
    double &newLoglk,
    bool blo,
    std::unordered_map<unsigned int, double> *treeHashScore = nullptr
    );


  /**
   *  Test all input SPR moves in parallel, and fill sortedBetterMoves
   *  with all the SPR moves that improve the likelihood, sorted by
   *  decreasing likelihood
   *
   *  @param jointTree The current tree
   *  @param allMoves The moves to test
   *  @param sortedBetterMoves The output better moves
   *  @param blo Do we apply branch length optimization?
   *  @param check Safety check for debugging purpose (it's expensive)
   *
   *  Parallelization: this function is parallelized over the moves 
   *  to test. 
   *
   */
  static bool findBetterMoves(JointTree &jointTree,
    std::vector<std::shared_ptr<SPRMove> > &allMoves,
    std::vector<std::shared_ptr<SPRMove> > &sortedBetterMoves,
    bool blo,
    bool check);


  /**
   *
   *  Find the best SPR move for a given prune index under a given radius
   *
   *
   *
   *  Parallelization: this function is local to one rank
   */
  static bool diggBestMoveFromPrune(JointTree &jointTree,
    std::unordered_map<unsigned int, double> &treeHashScores,
    unsigned int pruneIndex,
    unsigned int maxRadius,
    unsigned int additionalRadius,
    double &bestLL,
    double &bestLLAmongPrune,
    bool blo,
    SPRMove &bestMove);

};

