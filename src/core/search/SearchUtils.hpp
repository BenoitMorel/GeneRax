#pragma once

#include <search/Moves.hpp>

#include <memory>



class JointTree;

class SearchUtils {
public:
  static void testMove(JointTree &jointTree,
    std::shared_ptr<Move> move,
    double initialReconciliationLoglk,
    double initialLibpllLoglk,
    double &averageReconciliationDiff,
    double &newLoglk,
    bool blo,
    bool check
    );
 
  static bool findBestMove(JointTree &jointTree,
    std::vector<std::shared_ptr<Move> > &allMoves,
    double &bestLoglk,
    int &bestMoveIndex,
    bool blo,
    bool check);
};

