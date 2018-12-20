#include "SearchUtils.hpp"
#include <treeSearch/JointTree.hpp>
#include <Arguments.hpp>
#include <ParallelContext.hpp>
#include <Logger.hpp>



void SearchUtils::testMove(JointTree &jointTree,
    shared_ptr<Move> move,
    double initialReconciliationLoglk,
    double initialLibpllLoglk,
    double &averageReconciliationDiff,
    double &newLoglk
    )
{
  double initialLoglk = initialReconciliationLoglk + initialLibpllLoglk;
  jointTree.applyMove(move); 
  double recLoglk = jointTree.computeReconciliationLoglk();
  double improvement = recLoglk - initialReconciliationLoglk;
  averageReconciliationDiff *= 50;
  averageReconciliationDiff += improvement;
  averageReconciliationDiff /= 51;
  if (improvement < averageReconciliationDiff) {
    jointTree.rollbackLastMove();
    return;
  }
  jointTree.optimizeMove(move);
  newLoglk = recLoglk +  jointTree.computeLibpllLoglk(true);
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
  double initialReconciliationLoglk = jointTree.computeReconciliationLoglk();
  double initialLibpllLoglk = jointTree.computeLibpllLoglk();
  double averageReconciliationDiff = 0;
  assert(fabs(initialLoglk 
        - (initialReconciliationLoglk + initialLibpllLoglk)) < 0.000000001);
  int begin = ParallelContext::getBegin(allMoves.size());
  int end = ParallelContext::getEnd(allMoves.size());
  for (int i = begin; i < end; ++i) {
    auto move = allMoves[i];
    double loglk = bestLoglk;
    SearchUtils::testMove(jointTree, move, 
        initialReconciliationLoglk,
        initialLibpllLoglk, 
        averageReconciliationDiff,
        loglk);
    if (loglk > bestLoglk + 0.000000001) {
      bestLoglk = loglk;
      bestMoveIndex = i;
    }
  }
  int bestRank = 0;
  ParallelContext::getBestLL(bestLoglk, bestRank);
  ParallelContext::broadcoastInt(bestRank, bestMoveIndex);
  return bestMoveIndex != -1;
}

