#include "SPRSearch.h"
#include "Arguments.hpp"
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <ale/containers/GeneMap.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>

struct SPRMoveDesc {
  SPRMoveDesc(int prune, int regraft, const vector<int> &edges):
    pruneIndex(prune), regraftIndex(regraft), path(edges) {}
  int pruneIndex;
  int regraftIndex;
  vector<int> path;
};

void queryPruneIndicesRec(pll_unode_t * node,
                               vector<int> &buffer)
{
  if (node->next) {
    queryPruneIndicesRec(node->next->back, buffer);
    queryPruneIndicesRec(node->next->next->back, buffer);
    buffer.push_back(node->node_index);
  }
}

void getAllPruneIndices(JointTree &tree, vector<int> &allNodeIndices) {
  auto treeinfo = tree.getTreeInfo();
  vector<int> allNodes;
  queryPruneIndicesRec(treeinfo->root->back, allNodeIndices);
  queryPruneIndicesRec(treeinfo->root, allNodeIndices);
}


bool sprYeldsSameTree(pll_unode_t *n1, pll_unode_t *n2)
{
  return (n2 == n1) || (n2 == n1->next) || (n2 == n1->next->next)
    || (n2->back == n1) || (n2->back == n1->next) || (n2->back == n1->next->next);
}

bool isValidSPRMove(shared_ptr<pllmod_treeinfo_t> treeinfo, pll_unode_s *prune, pll_unode_s *regraft) {
  return !sprYeldsSameTree(prune, regraft);
}

bool testSPRMove(JointTree &jointTree,
    const SPRMoveDesc &move,
    double bestLoglk,
    double &newLoglk,
    shared_ptr<Move> &newMove)
{
  int pruneIndex = move.pruneIndex;
  int regraftIndex = move.regraftIndex;
  if (!isValidSPRMove(jointTree.getTreeInfo(), jointTree.getNode(pruneIndex), jointTree.getNode(regraftIndex))) {
    return false;
  }
  bool res = false;
  double initialLoglk = 0.0;
  if (Arguments::check) {
    initialLoglk = jointTree.computeJointLoglk();
  }
  newMove = Move::createSPRMove(pruneIndex, regraftIndex, move.path);
  jointTree.applyMove(newMove);
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
  return newLoglk > bestLoglk;
}

void getRegraftsRec(int pruneIndex, pll_unode_t *regraft, int maxRadius, vector<int> &path, vector<SPRMoveDesc> &moves)
{
  moves.push_back(SPRMoveDesc(pruneIndex, regraft->node_index, path));
  if (path.size() < maxRadius && regraft->next) {
    path.push_back(regraft->node_index);
    getRegraftsRec(pruneIndex, regraft->next->back, maxRadius, path, moves);
    getRegraftsRec(pruneIndex, regraft->next->next->back, maxRadius, path, moves);
    path.pop_back();
  }
}

void getRegrafts(JointTree &jointTree, int pruneIndex, int maxRadius, vector<SPRMoveDesc> &moves) 
{
  pll_unode_t *pruneNode = jointTree.getNode(pruneIndex);
  vector<int> path;
  getRegraftsRec(pruneIndex, pruneNode->next->back, maxRadius, path, moves);
  getRegraftsRec(pruneIndex, pruneNode->next->next->back, maxRadius, path, moves);
}

bool SPRSearch::applySPRRound(ParallelJointTree &jointTree, int radius, double &bestLoglk) {
  if (Arguments::verbose) {
    cout << "SPR ROUND WITH RADIUS " << radius << endl;
  }
    shared_ptr<Move> bestMove(0);
  vector<int> allNodes;
  getAllPruneIndices(jointTree.getThreadInstance(), allNodes);
  const size_t edgesNumber = jointTree.getThreadInstance().getTreeInfo()->tree->edge_count;
  bool foundBetterMove = false;
 
  vector<SPRMoveDesc> movesToExplore;
  for (int i = 0; i < allNodes.size(); ++i) {
      int pruneIndex = allNodes[i];
      getRegrafts(jointTree.getThreadInstance(), pruneIndex, radius, movesToExplore);
  }

  #pragma omp parallel for num_threads(jointTree.getThreadsNumber())
  for (int i = 0; i < movesToExplore.size(); ++i) {
    if (true) {
        double newLoglk;
        shared_ptr<Move> newMove;
        if (testSPRMove(jointTree.getThreadInstance(), movesToExplore[i], bestLoglk, newLoglk, newMove)) {
          #pragma omp critical
          if (bestLoglk < newLoglk) {
            foundBetterMove = true;
            bestMove = newMove;
            bestLoglk = newLoglk;
            if (Arguments::verbose) {
              cout << "found a better move with loglk " << newLoglk << endl;
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


void SPRSearch::applySPRSearch(ParallelJointTree &jointTree)
{ 
  jointTree.getThreadInstance().printLoglk();
  double startingLoglk = jointTree.getThreadInstance().computeJointLoglk();
  double bestLoglk = startingLoglk;
  while (applySPRRound(jointTree, 1, bestLoglk)) {}
  jointTree.optimizeParameters();
  while (applySPRRound(jointTree, 2, bestLoglk)) {}
  jointTree.optimizeParameters();
  while (applySPRRound(jointTree, 5, bestLoglk)) {}
}

