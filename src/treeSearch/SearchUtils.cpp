#include "SearchUtils.hpp"

#include <treeSearch/JointTree.h>
#include <Arguments.hpp>
#include <ParallelContext.hpp>
#include <Logger.hpp>

void SearchUtils::testMove(JointTree &jointTree,
    shared_ptr<Move> move,
    double initialLoglk,
    double &newLoglk
    )
{
  jointTree.applyMove(move);
  newLoglk = jointTree.computeJointLoglk();
  jointTree.rollbackLastMove();
  if(Arguments::check && fabs(initialLoglk - jointTree.computeJointLoglk()) > 0.000001) {
    cerr.precision(17);
    cerr << "rollback lead to different likelihoods: " << initialLoglk
      << " " << jointTree.computeJointLoglk() << endl;
    cerr << " rank " << ParallelContext::getRank() << endl;
    cerr << "Move: " << *move << endl;
    exit(1);
  }
}
  
bool SearchUtils::findBestMove(JointTree &jointTree,
    vector<shared_ptr<Move> > &allMoves,
    double &bestLoglk,
    int &bestMoveIndex)
{
  bestMoveIndex = -1;
  double initialLoglk = bestLoglk; //jointTree.computeJointLoglk();
  int begin = ParallelContext::getBegin(allMoves.size());
  int end = ParallelContext::getEnd(allMoves.size());
  for (int i = begin; i < end; ++i) {
    auto move = allMoves[i];
    double loglk;
    SearchUtils::testMove(jointTree, move, initialLoglk, loglk);
    if (loglk > bestLoglk) {
      bestLoglk = loglk;
      bestMoveIndex = i;
    }
  }
  int bestRank = 0;
  ParallelContext::getBestLL(bestLoglk, bestRank);
  ParallelContext::broadcoastInt(bestRank, bestMoveIndex);
  return bestMoveIndex != -1;
}
