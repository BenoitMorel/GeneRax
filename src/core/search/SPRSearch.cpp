#include <search/SPRSearch.hpp>
#include <trees/JointTree.hpp>
#include <search/Moves.hpp>
#include <search/SearchUtils.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>

#include <unordered_set>
#include <array>

struct SPRMoveDesc {
  SPRMoveDesc(unsigned int prune, unsigned int regraft, const std::vector<unsigned int> &edges):
    pruneIndex(prune), regraftIndex(regraft), path(edges) {}
  unsigned int pruneIndex;
  unsigned int regraftIndex;
  std::vector<unsigned int> path;
};

static void queryPruneIndicesRec(corax_unode_t * node,
                               std::vector<unsigned int> &buffer)
{
  assert(node);
  if (node->next) {
    queryPruneIndicesRec(node->next->back, buffer);
    queryPruneIndicesRec(node->next->next->back, buffer);
    buffer.push_back(node->node_index);
    if (node->back->next) {
      buffer.push_back(node->back->node_index);
    }
  }
}

static void getAllPruneIndices(JointTree &tree, std::vector<unsigned int> &allNodeIndices) {
  auto treeinfo = tree.getTreeInfo();
  for (unsigned int i = 0; i < treeinfo->subnode_count; ++i) {
    if (treeinfo->subnodes[i]->next) {
      allNodeIndices.push_back(treeinfo->subnodes[i]->node_index);
    }
  }
}


static bool sprYeldsSameTree(corax_unode_t *p, corax_unode_t *r)
{
  assert(p);
  assert(r);
  return (r == p) || (r == p->next) || (r == p->next->next)
    || (r == p->back) || (r == p->next->back) || (r == p->next->next->back);
}

static bool isValidSPRMove(corax_unode_s *prune, corax_unode_s *regraft) {
  assert(prune);
  assert(regraft);
  return !sprYeldsSameTree(prune, regraft);
}

static bool isSPRMoveValid(PLLUnrootedTree &tree,
    corax_unode_t *prune, 
    corax_unode_t *regraft)
{
  // regraft should not be a child of prune
  auto pruneChildren = tree.getPostOrderNodesFrom(prune->back);
  for (auto child: pruneChildren) {
    if (regraft == child) {
      return false;
    }
    if (regraft->back == child) {
      return false;
    }
  }
  return !sprYeldsSameTree(prune, regraft);
}



static void getRegraftsRec(unsigned int pruneIndex, 
    corax_unode_t *regraft, 
    int maxRadius, 
    double supportThreshold, 
    std::vector<unsigned int> &path, 
    std::vector<SPRMoveDesc> &moves)
{
  assert(regraft);
  double bootstrapValue = (nullptr == regraft->label) ? 0.0 : std::atof(regraft->label);
  if (supportThreshold >= 0.0 && bootstrapValue > supportThreshold) {
    return;
  }
  if (path.size()) {
    moves.push_back(SPRMoveDesc(pruneIndex, regraft->node_index, path));
  }
  if (static_cast<int>(path.size()) < maxRadius && regraft->next) {
    path.push_back(regraft->node_index);
    getRegraftsRec(pruneIndex, regraft->next->back, maxRadius, supportThreshold, path, moves);
    getRegraftsRec(pruneIndex, regraft->next->next->back, maxRadius, supportThreshold, path, moves);
    path.pop_back();
  }
}

static void getRegrafts(JointTree &jointTree, unsigned int pruneIndex, int maxRadius, std::vector<SPRMoveDesc> &moves) 
{
  corax_unode_t *pruneNode = jointTree.getNode(pruneIndex);
  std::vector<unsigned int> path;
  auto supportThreshold = jointTree.getSupportThreshold();
  getRegraftsRec(pruneIndex, pruneNode->next->back, maxRadius, supportThreshold, path, moves);
  getRegraftsRec(pruneIndex, pruneNode->next->next->back, maxRadius, supportThreshold, path, moves);
}

static void addInvolvedNode(corax_unode_t *node, 
    std::unordered_set<corax_unode_t *> &involved)
{
  involved.insert(node);
  if (node->next) {
    involved.insert(node->next);
    involved.insert(node->next->next);
  }
}

static bool wasInvolved(corax_unode_t *node,
    std::unordered_set<corax_unode_t *> &involved)
{
  return involved.find(node) != involved.end();
}

bool SPRSearch::applySPRRound(JointTree &jointTree, int radius, double &bestLoglk, bool blo) {
  std::vector<unsigned int> allNodes;
  getAllPruneIndices(jointTree, allNodes);
  std::vector<SPRMoveDesc> potentialMoves;
  std::vector<std::shared_ptr<SPRMove> > allMoves;
  for (unsigned int i = 0; i < allNodes.size(); ++i) {
      auto pruneIndex = allNodes[i];
      getRegrafts(jointTree, pruneIndex, radius, potentialMoves);
  }
  std::vector<std::array<bool, 2> > redundantNNIMoves(
      jointTree.getTreeInfo()->subnode_count, 
      std::array<bool, 2>{{false,false}});
  for (auto &move: potentialMoves) {
    auto pruneIndex = move.pruneIndex;
    auto regraftIndex = move.regraftIndex;
    if (!isValidSPRMove(jointTree.getNode(pruneIndex), jointTree.getNode(regraftIndex))) {
      continue;
    }
    // in case of a radius 1, SPR moves are actually NNI moves
    // the traversal algorithm produces redundant moves, so we 
    // get rid of them here
    if (move.path.size() == 1) { 
      auto nniEdge = jointTree.getNode(move.path[0]);
      bool isPruneNext = nniEdge->back->next->node_index == pruneIndex;
      bool isRegraftNext = nniEdge->next->back->node_index == regraftIndex;
      auto nniType = static_cast<unsigned int>(isPruneNext == isRegraftNext);
      auto nniBranchIndex = std::min(nniEdge->node_index, nniEdge->back->node_index);
      if (redundantNNIMoves[nniBranchIndex][nniType]) {
        continue;
      }
      redundantNNIMoves[nniBranchIndex][nniType] = true; 
    }
    allMoves.push_back(SPRMove::createSPRMove(pruneIndex, regraftIndex, move.path));
  }
  Logger::info << "Start SPR round " 
    << "(std::hash=" << jointTree.getUnrootedTreeHash() << ", (best ll=" 
    << bestLoglk << ", radius=" << radius << ", possible moves: " << allMoves.size() << ")"
    << std::endl;
  std::vector<std::shared_ptr<SPRMove> > betterMoves;
  double initialLL = bestLoglk; 
  SearchUtils::findBetterMoves(jointTree, 
      allMoves,
      betterMoves,
      blo,
      jointTree.isSafeMode());
  bool foundBetterMove = false;
  std::unordered_set<corax_unode_t *> involved;
  //for (auto &move: betterMoves) {
  double currentLL = jointTree.computeJointLoglk();
  for (unsigned int i = 0; i < betterMoves.size(); ++i) {
    auto move = betterMoves[i];
    auto prune = jointTree.getNode(move->getPruneIndex());
    auto regraft = jointTree.getNode(move->getRegraftIndex());
    assert(ParallelContext::isIntEqual(i));
    assert(ParallelContext::isIntEqual(move->getPruneIndex()));
    if (wasInvolved(prune, involved) || wasInvolved(regraft, involved)) {
      Logger::info << "involved, abort" << std::endl;
      continue;
    }
    if (!isSPRMoveValid(jointTree.getGeneTree(), prune, regraft)) {
      Logger::info << "invalid, abort" << std::endl;
      continue;
    }
    addInvolvedNode(prune, involved);
    addInvolvedNode(regraft, involved);
    move->updatePath(jointTree);
    jointTree.applyMove(*move);
    if (blo) {
      jointTree.optimizeMove(*move);
    }
    double ll = jointTree.computeJointLoglk();
    if (ll < currentLL) {
      Logger::info << "ROLLBACK" << std::endl;
      jointTree.rollbackLastMove();
      ll = jointTree.computeJointLoglk();
    } else {
      foundBetterMove = true;
      currentLL = ll;
      bestLoglk = ll;
    }
    Logger::info << move->getScore() << " " << ll << std::endl;
  }
  double epsilon = fabs(bestLoglk) > 10000.0 ? 0.5 : 0.001; 
  return foundBetterMove && (bestLoglk - initialLL > epsilon);
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
  jointTree.optimizeParameters(true, false);
  bestLoglk = jointTree.computeJointLoglk();
  while (applySPRRound(jointTree, 5, bestLoglk)) {}
}

