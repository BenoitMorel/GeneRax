#pragma once

#include <search/Moves.hpp>

#include <memory>



class JointTree;

class SearchUtils {
public:
  static void testMove(JointTree &jointTree,
    SPRMove &move,
    double initialReconciliationLoglk,
    double initialLibpllLoglk,
    double &averageReconciliationDiff,
    double &newLoglk,
    bool blo,
    bool check
    );
 
  static bool findBestMove(JointTree &jointTree,
    std::vector<std::shared_ptr<SPRMove> > &allMoves,
    double &bestLoglk,
    unsigned int &bestMoveIndex,
    bool blo,
    bool check);
  
  static bool findBetterMoves(JointTree &jointTree,
    std::vector<std::shared_ptr<SPRMove> > &allMoves,
    std::vector<std::shared_ptr<SPRMove> > &sortedBetterMoves,
    bool blo,
    bool check);

};

