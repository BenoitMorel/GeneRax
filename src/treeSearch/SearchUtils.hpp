#pragma once

#include <treeSearch/Moves.hpp>

#include <memory>

using namespace std;

class JointTree;

class SearchUtils {
public:
  static void testMove(JointTree &jointTree,
    shared_ptr<Move> move,
    double initialReconciliationLoglk,
    double initialLibpllLoglk,
    double &averageReconciliationDiff,
    double &newLoglk
    );
 
  static bool findBestMove(JointTree &jointTree,
    vector<shared_ptr<Move> > &allMoves,
    double &bestLoglk,
    int &bestMoveIndex);
};

