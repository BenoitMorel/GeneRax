#include "NNISearch.h"
#include <algorithm>
#include <omp.h>
bool testNNIMove(JointTree &jointTree,
    int nodeIndex,
    int nniType,
    double bestLoglk,
    double &newLoglk,
    shared_ptr<Move> &newMove
    )
{
  //double initialLoglk = jointTree.computeJointLoglk();
  newMove = Move::createNNIMove(nodeIndex, nniType, true);
  jointTree.applyMove(newMove);
  newLoglk = jointTree.computeLibpllLoglk();
  if (newLoglk > bestLoglk) {
    newLoglk += jointTree.computeALELoglk();
  }
  jointTree.rollbackLastMove();
  //assert(initialLoglk == jointTree.computeJointLoglk());
  return newLoglk > bestLoglk;
}

bool NNISearch::applyNNIRound(AbstractJointTree &jointTree, double &bestLoglk) {
  cout << "Start NNI Round" << endl;
  shared_ptr<Move> bestMove(0);
  vector<int> allNodes;
  jointTree.getThreadInstance().getAllNodeIndices(allNodes);
  bool foundBetterMove = false;
  #pragma omp parallel for num_threads(jointTree.getThreadsNumber())
  for (int j = 0; j < allNodes.size(); ++j) {
    if (!foundBetterMove) { 
      for (int moveType = 0; moveType < 2; ++moveType) {
        shared_ptr<Move> newMove;
        double loglk;
        if (testNNIMove(jointTree.getThreadInstance(), allNodes[j], moveType, bestLoglk, loglk, newMove)) {
          #pragma omp critical
          if (!foundBetterMove) {
            foundBetterMove = true;
            bestMove = newMove;
            bestLoglk = loglk;
            cout << "found a better move with loglk " << loglk << endl;
          }
        }
      }
    }
  }
  if (foundBetterMove) {
    jointTree.applyMove(bestMove);
  }
  return foundBetterMove;
}


void NNISearch::applyNNISearch(AbstractJointTree &jointTree)
{
  jointTree.getThreadInstance().printLoglk();
  double bestLoglk = jointTree.getThreadInstance().computeJointLoglk();
  bool foundBetterMove = true;
  while (applyNNIRound(jointTree, bestLoglk)) {
  }
}
