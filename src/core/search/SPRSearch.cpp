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

static void getAllPotentialMoves(JointTree &jointTree,
    int radius,
    std::vector<std::shared_ptr<SPRMove> > &allMoves)
{
  std::vector<unsigned int> allNodes;
  getAllPruneIndices(jointTree, allNodes);
  std::vector<SPRMoveDesc> potentialMoves;
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

}

static bool sequentiallyApplyBetterMoves(JointTree &jointTree, 
    const std::vector<std::shared_ptr<SPRMove> > &betterMoves,
    bool blo,
    double &bestLoglk)
{
  Logger::timed << "GeneRax will now test and apply all the " << betterMoves.size() << " potential good moves one after each other" << std::endl;
  bool foundBetterMove = false;
  std::unordered_set<corax_unode_t *> involved;
  for (auto move: betterMoves) {
    auto prune = jointTree.getNode(move->getPruneIndex());
    auto regraft = jointTree.getNode(move->getRegraftIndex());
    if (wasInvolved(prune, involved) || wasInvolved(regraft, involved)) {
      continue;
    }
    if (!isSPRMoveValid(jointTree.getGeneTree(), prune, regraft)) {
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
    if (ll < bestLoglk) {
      Logger::info << "Move rejected" << std::endl;
      jointTree.rollbackLastMove();
      ll = jointTree.computeJointLoglk();
    } else {
      foundBetterMove = true;
      bestLoglk = ll;
      Logger::timed << "\tApplying move, ll = " << ll << std::endl;
    }
  }
  return foundBetterMove;
}

static bool simultaneouslyApplyBetterMoves(JointTree &jointTree, 
    const std::vector<std::shared_ptr<SPRMove> > &betterMoves,
    bool blo,
    double &bestLoglk,
    std::vector<std::shared_ptr<SPRMove> > &movesToRetry)
{
  Logger::timed << "GeneRax will now test and apply all non-conflicting " << betterMoves.size() << " potential good moves simultaneously" << std::endl;
  std::unordered_set<corax_unode_t *> involved;
  std::vector<std::shared_ptr<SPRMove> >  appliedMoves;
  for (auto move: betterMoves) {
    auto prune = jointTree.getNode(move->getPruneIndex());
    auto regraft = jointTree.getNode(move->getRegraftIndex());
    if (!isSPRMoveValid(jointTree.getGeneTree(), prune, regraft)) {
      continue;
    }
    auto p = prune;
    auto r = regraft;
    std::vector<corax_unode_t *> branches;
    PLLUnrootedTree::orientTowardEachOther(&p, &r, branches);
    bool skip = false;
    for (auto b: branches) {
      if (wasInvolved(b, involved)) {
        skip = true;
      }
    }
    if (skip) {
      movesToRetry.push_back(move);
      continue;
    }
    for (auto b: branches) {
      addInvolvedNode(b, involved);
      addInvolvedNode(b->back, involved);
    }
    move->updatePath(jointTree);
    jointTree.applyMove(*move);
    appliedMoves.push_back(move);
    if (blo) {
      jointTree.reOptimizeMove(*move);
    }
  }
  double initialLoglk = bestLoglk;
  bestLoglk = jointTree.computeJointLoglk();
  Logger::timed << "After simultaneously applying the moves, ll = " << bestLoglk << std::endl;
  if (bestLoglk < initialLoglk) {
    Logger::info <<"Rollback all moves!" << std::endl;
    for (auto move: appliedMoves) {
      jointTree.rollbackLastMove();
      movesToRetry.push_back(move);
    }
    bestLoglk = jointTree.computeJointLoglk();
    return false;
  }
  return true;
}


bool SPRSearch::applySPRRound(JointTree &jointTree, int radius, double &bestLoglk, bool blo) {
  std::vector<std::shared_ptr<SPRMove> > allMoves;
  getAllPotentialMoves(jointTree, radius, allMoves);
  Logger::timed << "Start SPR round " 
    << "(std::hash=" << jointTree.getUnrootedTreeHash() << ", (best ll=" 
    << bestLoglk << ", radius=" << radius << ", possible moves: " << allMoves.size() << ")"
    << std::endl;
  std::vector<std::shared_ptr<SPRMove> > betterMoves;
  Logger::info << "Searching all potential good moves..." << std::endl;
  SearchUtils::findBetterMoves(jointTree, 
      allMoves,
      betterMoves,
      blo,
      jointTree.isSafeMode());
  bestLoglk = jointTree.computeJointLoglk();
  double initialLL = bestLoglk; 
  bool foundBetterMove = false;
  
  if (betterMoves.size() > 10) {
    std::vector<std::shared_ptr<SPRMove> > movesToRetry;
    foundBetterMove |= simultaneouslyApplyBetterMoves(jointTree, 
        betterMoves,
        blo,
        bestLoglk,
        movesToRetry);
    if (movesToRetry.size()) {
      Logger::info << "Now sequentially applying the " << movesToRetry.size() << " moves that could not be applied simultaneously..." << std::endl;
      foundBetterMove |= sequentiallyApplyBetterMoves(jointTree, 
          movesToRetry,
          blo,
          bestLoglk);
    }
  } else {
    foundBetterMove |= sequentiallyApplyBetterMoves(jointTree, 
        betterMoves,
        blo,
        bestLoglk);
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

