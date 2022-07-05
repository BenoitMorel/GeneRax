#pragma once

#include <search/Moves.hpp>
#include <maths/AverageStream.hpp>
#include <memory>

#include <unordered_map>

class JointTree;


/**
 *  Store information about the previously tested move likelihoods 
 *  Provide a stopping criteria for the SPR search, to stop before
 *  we reach the maximum SPR radius if the current likelihood is too
 *  low
 */
struct DiggStopper {
  DiggStopper(): worseDiffR1(0.0), streamR1(50), stream(50), ok(0), ko(0)  {}
  double worseDiffR1;
  AverageStream streamR1;
  AverageStream stream;
  unsigned int ok;
  unsigned int ko;

  void treatDiff(double diff, unsigned int radius) {
    if (diff > 0.1) {
      return;
    }
    if (radius == 1) {
      worseDiffR1 = std::min(diff, worseDiffR1);
      streamR1.addValue(diff);
    }
    stream.addValue(diff);
  }

  bool doStop(double diff, unsigned int) {
    //return false;
    bool res = stream.isSignificant() && (diff < streamR1.getAverage());
    //bool res = stream.isSignificant() && (diff < worseDiffR1 / 2.0);
    if (res) {
      ko += 1;
    } else {
      ok += 1;
    }
    return res;
  }
};


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
    bool blo);

  /**
   *
   *  Find the best SPR move for a given prune index under a given radius
   *  
   *  @param jointTree The current tree
   *  @param treeHashScore A buffer storing likelihoods to avoid
   *     computing the likelihood of the same tree twice
   *     This object should be cleared when the DTL rates are
   *     optimized
   *  @param stopper Criterion to stop the search earlier
   *  @param pruneIndex The node index of the node to prune
   *  @param maxRadius Maximum SPR radius
   *  @param additionalRadius Value added to maxRadius when a better
   *    move is found to search further
   *  @param bestLL the best likelihood (should be initialized
   *    by the caller)
   *  @param bestLLAmongPrune the best likelihood among all 
   *    tested moves (can be worse than the initial 
   *    tree likelihood). Initial value will be ignored
   *  @param blo Should we apply branch length opt
   *  @param bestMove Best move (even if it has a lower 
   *    likelihood than the initial tree)
   *
   *  Parallelization: this function is local to one rank
   */
  static bool diggBestMoveFromPrune(JointTree &jointTree,
    std::unordered_map<unsigned int, double> &treeHashScores,
    DiggStopper &stopper,
    unsigned int pruneIndex,
    unsigned int maxRadius,
    unsigned int additionalRadius,
    double &bestLL,
    double &bestLLAmongPrune,
    bool blo,
    SPRMove &bestMove);

};

