#pragma once

#include <search/Moves.hpp>

#include <memory>



class JointTree;

class SearchUtils {
public:
  static void testMove(JointTree &jointTree,
    Move &move,
    double initialReconciliationLoglk,
    double initialLibpllLoglk,
    double &averageReconciliationDiff,
    double &newLoglk,
    bool blo,
    bool check
    );
 
  static bool findBestMove(JointTree &jointTree,
    std::vector<std::unique_ptr<Move> > &allMoves,
    double &bestLoglk,
    unsigned int &bestMoveIndex,
    bool blo,
    bool check);
};

