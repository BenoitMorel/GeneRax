#include "NNISearch.h"
#include <algorithm>
#include "Arguments.hpp"
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <ale/containers/GeneMap.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>
#include <ale/tools/SpeciesGeneMapper.h>
#include <treeSearch/SearchUtils.hpp>
#include <Logger.hpp>

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


bool NNISearch::applyNNIRound(JointTree &jointTree, double &bestLoglk, int bloRadius) {
  shared_ptr<Move> bestMove(0);
  vector<int> allNodes;
  vector<shared_ptr<Move> > allMoves;
  getAllNNIIndices(jointTree, allNodes);
  for (int j = 0; j < allNodes.size(); ++j) {
      for (int moveType = 0; moveType < 2; ++moveType) {
        allMoves.push_back(Move::createNNIMove(allNodes[j], moveType, true, bloRadius));
      }
  }
  Logger::timed << "Start NNI Round (best ll=" << bestLoglk << ", " << allMoves.size() << " moves to try)" << endl;
  int bestMoveIndex = -1;
  bool foundBetterMove = SearchUtils::findBestMove(jointTree, allMoves, bestLoglk, bestMoveIndex); 
  if (foundBetterMove) {
    jointTree.applyMove(allMoves[bestMoveIndex]);
  }
  return foundBetterMove;
}


void NNISearch::applyNNISearch(JointTree &jointTree)
{
  jointTree.printLoglk();
  double bestLoglk = jointTree.computeJointLoglk();
  bool foundBetterMove = true;
  while (applyNNIRound(jointTree, bestLoglk, 1)) {}
  jointTree.optimizeParameters();
  bestLoglk = jointTree.computeJointLoglk();
  while (applyNNIRound(jointTree, bestLoglk, 3)) {}
}
