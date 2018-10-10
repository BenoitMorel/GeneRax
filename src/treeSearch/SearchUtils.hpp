#ifndef JOINT_TREE_SEARCH_SEARCH_UTILS_HPP_
#define JOINT_TREE_SEARCH_SEARCH_UTILS_HPP_

#include <memory>
#include <treeSearch/Moves.h>

using namespace std;

class JointTree;

class SearchUtils {
public:
  static void testMove(JointTree &jointTree,
    shared_ptr<Move> move,
    double bestLoglk,
    double &newLoglk
    );
 
  static bool findBestMove(JointTree &jointTree,
    vector<shared_ptr<Move> > &allMoves,
    double &bestLoglk,
    int &bestMoveIndex);

};

#endif
