#include "NNISearch.h"
#include <algorithm>
#include <omp.h>

void queryNNIIndicesRec(pll_unode_t * node,
                               vector<int> &buffer,
                               bool include_Node = true)
{
  if (node->next) {
    queryNNIIndicesRec(node->next->back, buffer);
    queryNNIIndicesRec(node->next->next->back, buffer);
    if (include_Node && node->back->next) {
        buffer.push_back(node->node_index);
    }
  }
}

void getAllNNIIndices(JointTree &tree, vector<int> &allNodeIndices) {
  auto treeinfo = tree.getTreeInfo();
  vector<int> allNodes;
  queryNNIIndicesRec(treeinfo->root->back, allNodeIndices, false);
  queryNNIIndicesRec(treeinfo->root, allNodeIndices);
}


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
  getAllNNIIndices(jointTree.getThreadInstance(), allNodes);
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
