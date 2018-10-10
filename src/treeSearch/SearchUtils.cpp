#include "SearchUtils.hpp"

#include <treeSearch/JointTree.h>
#include <Arguments.hpp>

void SearchUtils::testMove(JointTree &jointTree,
    shared_ptr<Move> move,
    double bestLoglk,
    double &newLoglk
    )
{
  double initialLoglk = 0.0;
  if (Arguments::check) {
    jointTree.computeJointLoglk();
  }
  jointTree.applyMove(move);
  newLoglk = jointTree.computeLibpllLoglk();
  if (newLoglk > bestLoglk) {
    newLoglk += jointTree.computeALELoglk();
  }
  jointTree.rollbackLastMove();
  if(Arguments::check && fabs(initialLoglk - jointTree.computeJointLoglk()) > 0.000001) {
    cerr.precision(17);
    cerr << "rollback lead to different likelihoods: " << initialLoglk
      << " " << jointTree.computeJointLoglk() << endl;
    exit(1);
  }
}
  
bool SearchUtils::findBestMove(JointTree &jointTree,
    vector<shared_ptr<Move> > &allMoves,
    double &bestLoglk,
    int &bestMoveIndex)
{
  bestMoveIndex = -1;
  bool foundBetterMove = false;
  for (int i = 0; i < allMoves.size(); ++i) {
    auto move = allMoves[i];
    if (!foundBetterMove) {
      double loglk;
      SearchUtils::testMove(jointTree, move, bestLoglk, loglk);
      if (loglk > bestLoglk) {
        if (!foundBetterMove) {
          foundBetterMove = true;
          bestLoglk = loglk;
          bestMoveIndex = i;
          if (Arguments::verbose) {
            cout << "found a better move with loglk " << loglk << endl;
          }
        }
      }
    }
  }
  return foundBetterMove;
}

