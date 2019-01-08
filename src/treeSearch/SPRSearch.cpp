#include <treeSearch/SPRSearch.hpp>
#include <treeSearch/JointTree.hpp>
#include <treeSearch/Moves.hpp>
#include <treeSearch/SearchUtils.hpp>
#include <Logger.hpp>
#include <Arguments.hpp>
#include <ParallelContext.hpp>

#include <unordered_set>
#include <array>

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
    if (node->back->next) {
      buffer.push_back(node->back->node_index);
    }
  }
}

void getAllPruneIndices(JointTree &tree, vector<int> &allNodeIndices) {
  auto treeinfo = tree.getTreeInfo();
  for (int i = 0; i < treeinfo->subnode_count; ++i) {
    if (treeinfo->subnodes[i]->next) {
      allNodeIndices.push_back(treeinfo->subnodes[i]->node_index);
    }
  }
}


bool sprYeldsSameTree(pll_unode_t *p, pll_unode_t *r)
{
  return (r == p) || (r == p->next) || (r == p->next->next)
    || (r == p->back) || (r == p->next->back) || (r == p->next->next->back);
}

bool isValidSPRMove(shared_ptr<pllmod_treeinfo_t> treeinfo, pll_unode_s *prune, pll_unode_s *regraft) {
  return !sprYeldsSameTree(prune, regraft);
}



void getRegraftsRec(int pruneIndex, pll_unode_t *regraft, int maxRadius, vector<int> &path, vector<SPRMoveDesc> &moves)
{
  if (path.size()) {
    moves.push_back(SPRMoveDesc(pruneIndex, regraft->node_index, path));
  }
  if (path.size() < maxRadius && regraft->next) {
    path.push_back(regraft->node_index);
    getRegraftsRec(pruneIndex, regraft->next->back, maxRadius, path, moves);
    getRegraftsRec(pruneIndex, regraft->next->next->back, maxRadius, path, moves);
    path.pop_back();
  }
}


void printPossibleMoves(JointTree &jointTree, vector<shared_ptr<Move> > &allMoves)
{
  //Logger::info << "Nodes: " << endl;
  //jointTree.printAllNodes(cout);
  Logger::info << "Possible moves from " << jointTree.getUnrootedTreeHash() << endl;
  unordered_set<int> hashs;
  for (auto move: allMoves) {
    jointTree.applyMove(move);
    auto hash = jointTree.getUnrootedTreeHash();
    hashs.insert(hash);
    jointTree.rollbackLastMove();
  }
  Logger::info << "Unique moves" << hashs.size() << endl;
}

void getRegrafts(JointTree &jointTree, int pruneIndex, int maxRadius, vector<SPRMoveDesc> &moves) 
{
  pll_unode_t *pruneNode = jointTree.getNode(pruneIndex);
  vector<int> path;
  getRegraftsRec(pruneIndex, pruneNode->next->back, maxRadius, path, moves);
  getRegraftsRec(pruneIndex, pruneNode->next->next->back, maxRadius, path, moves);
}

bool SPRSearch::applySPRRound(JointTree &jointTree, int radius, double &bestLoglk) {
  shared_ptr<Move> bestMove(0);
  vector<int> allNodes;
  getAllPruneIndices(jointTree, allNodes);
  const size_t edgesNumber = jointTree.getTreeInfo()->tree->edge_count;
 
  vector<SPRMoveDesc> potentialMoves;
  vector<shared_ptr<Move> > allMoves;
  for (int i = 0; i < allNodes.size(); ++i) {
      int pruneIndex = allNodes[i];
      getRegrafts(jointTree, pruneIndex, radius, potentialMoves);
  }
  vector<array<bool, 2> > redundantNNIMoves(jointTree.getTreeInfo()->subnode_count, array<bool, 2>{{false,false}});
  for (auto &move: potentialMoves) {
    int pruneIndex = move.pruneIndex;
    int regraftIndex = move.regraftIndex;
    if (!isValidSPRMove(jointTree.getTreeInfo(), jointTree.getNode(pruneIndex), jointTree.getNode(regraftIndex))) {
      continue;
    }
    // in case of a radius 1, SPR moves are actually NNI moves
    // the traversal algorithm produces redundant moves, so we 
    // get rid of them here
    if (move.path.size() == 1) { 
      auto nniEdge = jointTree.getNode(move.path[0]);
      bool isPruneNext = nniEdge->back->next->node_index == pruneIndex;
      bool isRegraftNext = nniEdge->next->back->node_index == regraftIndex;
      int nniType = (isPruneNext == isRegraftNext);
      int nniBranchIndex = min(nniEdge->node_index, nniEdge->back->node_index);
      if (redundantNNIMoves[nniBranchIndex][nniType]) {
        continue;
      }
      redundantNNIMoves[nniBranchIndex][nniType] = true; 
    }

    allMoves.push_back(Move::createSPRMove(pruneIndex, regraftIndex, move.path));
  }
  //printPossibleMoves(jointTree, allMoves);
  
  Logger::timed << "Start SPR round " 
    << "(hash=" << jointTree.getUnrootedTreeHash() << ", (best ll=" << bestLoglk << ", radius=" << radius << ", possible moves: " << allMoves.size() << ")"
    << endl;
  int bestMoveIndex = -1;
  bool foundBetterMove = SearchUtils::findBestMove(jointTree, allMoves, bestLoglk, bestMoveIndex); 
  if (foundBetterMove) {
    jointTree.applyMove(allMoves[bestMoveIndex]);
    jointTree.optimizeMove(allMoves[bestMoveIndex]);
    double ll = jointTree.computeJointLoglk();
    if (fabs(ll - bestLoglk) > 0.00000001) {
      cerr << ll << " " << bestLoglk << " " 
        << jointTree.computeJointLoglk() << endl;
    }
    assert(fabs(ll - bestLoglk) < 0.00000001);
  }
  return foundBetterMove;
}


void SPRSearch::applySPRSearch(JointTree &jointTree)
{ 
  jointTree.printLoglk();
  double startingLoglk = jointTree.computeJointLoglk();
  double bestLoglk = startingLoglk;

  while (applySPRRound(jointTree, 1, bestLoglk)) {}
  jointTree.optimizeParameters();
  bestLoglk = jointTree.computeJointLoglk();
  while (applySPRRound(jointTree, 1, bestLoglk)) {}
  jointTree.optimizeParameters(true, false);
  bestLoglk = jointTree.computeJointLoglk();
  while (applySPRRound(jointTree, 2, bestLoglk)) {}
  jointTree.optimizeParameters(true, false);
  bestLoglk = jointTree.computeJointLoglk();
  while (applySPRRound(jointTree, 3, bestLoglk)) {}
  //while (applySPRRound(jointTree, 50, bestLoglk)) {}
}

