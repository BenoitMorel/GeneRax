#include <search/SPRSearch.hpp>
#include <trees/JointTree.hpp>
#include <search/Moves.hpp>
#include <search/SearchUtils.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <unordered_map>
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
  /*
  auto rng = std::default_random_engine {};
  std::shuffle(std::begin(allNodeIndices), std::end(allNodeIndices), rng);
  */
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
  Logger::timed << "\tGeneRax will now test and apply all the " << betterMoves.size() << " potential good moves one after each other" << std::endl;
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
      Logger::info << "\tWorse likelihood, Move rejected" << std::endl;
      jointTree.rollbackLastMove();
      ll = jointTree.computeJointLoglk();
    } else {
      foundBetterMove = true;
      bestLoglk = ll;
      Logger::info << "\tApplying move, ll = " << ll << ", " << move->getPruneIndex() << " " << move->getRegraftIndex() << std::endl;
    }
  }
  return foundBetterMove;
}

static bool simultaneouslyApplyBetterMoves(JointTree &jointTree, 
    const std::vector<std::shared_ptr<SPRMove> > &betterMoves,
    bool blo,
    double &bestLoglk,
    std::unordered_set<corax_unode_t *> &involved,
    std::vector<std::shared_ptr<SPRMove> > &movesToRetry)
{
  double initialLoglk = bestLoglk;
  Logger::timed << "\tApplying a chunk of " << betterMoves.size() << " potential good moves simultaneously..." << std::endl;
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
  bestLoglk = jointTree.computeJointLoglk();
  if (bestLoglk < initialLoglk) {
    Logger::timed <<"\tWorse likelihood (" << bestLoglk << ")! Rollbacking all moves! " << bestLoglk << std::endl;
    for (auto move: appliedMoves) {
      jointTree.rollbackLastMove();
      movesToRetry.push_back(move);
    }
    bestLoglk = jointTree.computeJointLoglk();
    return false;
  }
  Logger::timed << "\tAfter applying the chunk of moves, ll = " << bestLoglk << std::endl;
  return true;
}

struct ScoredPrune {
  ScoredPrune(unsigned int i, double s): index(i), score(s) {}
  unsigned int index;
  double score;
};


struct less_than_prune
{
  inline bool operator() (const ScoredPrune &m1, 
      const ScoredPrune &m2)
  {
    if (m1.score != m2.score) {
      return (m1.score < m2.score);
    } 
    return m1.index < m2.index;
  }
};


//#define SIMPLE
static bool applyTheBetterMoves(JointTree &jointTree,
    const std::vector<std::shared_ptr<SPRMove> > &betterMoves,
    bool blo,
    double &bestLoglk)
{
  const unsigned int chunkSize = 10;
  const unsigned int maxSequential = 100;
  bool foundBetterMove = false;
  double initialLL = bestLoglk;
  Logger::timed << "Found " << betterMoves.size() << " potential better moves" << std::endl;
  if (betterMoves.size() > maxSequential) {
    Logger::timed << "GeneRax will try to apply them chunk by chunk..." << std::endl;
    unsigned int offset = 0;
    std::vector<std::shared_ptr<SPRMove> > movesToRetry;
    std::unordered_set<corax_unode_t *> involved;
    while (offset <= betterMoves.size()) {
      auto end = std::min<unsigned int>(offset + chunkSize, betterMoves.size());
      std::vector<std::shared_ptr<SPRMove> > chunk(betterMoves.begin() + offset,
          betterMoves.begin() + end);
      bool improved = simultaneouslyApplyBetterMoves(jointTree, 
        chunk,
        blo,
        bestLoglk,
        involved,
        movesToRetry);
      offset += chunkSize;
      if (improved) {
        // apply the moves that could not be applied simultaneously
        improved |= sequentiallyApplyBetterMoves(jointTree, 
            movesToRetry,
            blo,
            bestLoglk);
      } else {
        // apply all the moves one after each other
        improved |= sequentiallyApplyBetterMoves(jointTree, 
            chunk,
            blo,
            bestLoglk);
      }
      movesToRetry.clear();
      foundBetterMove |= improved;
      if (!improved) {
        break;
      }
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
  
static void synchronizeMoves(JointTree &jointTree, std::vector<std::shared_ptr<SPRMove> > &moves)
{
  std::vector<unsigned int> localPrune;    
  std::vector<unsigned int> localRegraft; 
  for (auto move: moves) {
    localPrune.push_back(move->getPruneIndex());
    localRegraft.push_back(move->getRegraftIndex());
  }
  std::vector<unsigned int> prune;    
  std::vector<unsigned int> regraft; 
  ParallelContext::concatenateHetherogeneousUIntVectors(localPrune, prune);
  ParallelContext::concatenateHetherogeneousUIntVectors(localRegraft, regraft);
  auto temp = moves;
  moves.clear();
  assert(prune.size() == regraft.size());
  std::vector<unsigned int> path;
  for (unsigned int i = 0; i < prune.size(); ++i) {
    moves.push_back(std::make_shared<SPRMove>(prune[i], regraft[i], path));
    moves.back()->updatePath(jointTree);
  }
}

bool SPRSearch::applySPRRoundDigg(JointTree &jointTree, int radius, bool blo) 
{
  std::vector<unsigned int> pruneIndices;
  getAllPruneIndices(jointTree, pruneIndices);
  unsigned int additionalRadius = 0;
  bool improved = false;
  double bestLL = jointTree.computeJointLoglk();
  std::unordered_map<unsigned int, double> treeHashScores;
  DiggStopper stopper;
  std::vector<std::shared_ptr<SPRMove> > betterMoves;
  std::vector<ScoredPrune> scoredPrunes;
  Logger::timed << "SPR Search with radius " << radius << ": trying " << pruneIndices.size() << " prune nodes" << std::endl;
  auto begin = ParallelContext::getBegin(pruneIndices.size());
  auto end = ParallelContext::getEnd(pruneIndices.size());
  for (unsigned int i = begin; i < end; ++i) {
    auto pruneIndex = pruneIndices[i];
    std::vector<unsigned int> path;
    SPRMove move(0, 0, path);
    double bestLLAmongPrune;
    double tempLL = bestLL;
    bool isBetter = SearchUtils::diggBestMoveFromPrune(jointTree, 
        treeHashScores,
        stopper,
        pruneIndex, 
        radius + 1,
        additionalRadius,
        tempLL,
        bestLLAmongPrune,
        blo,
        move);  
    //scoredPrunes.push_back(ScoredPrune(pruneIndex, bestLLAmongPrune - bestLL));
    if (isBetter) {
      betterMoves.push_back(std::make_shared<SPRMove>(move));
    }
  }
  synchronizeMoves(jointTree, betterMoves);
  improved |= applyTheBetterMoves(jointTree, 
    betterMoves,
    blo,
    bestLL);
  return improved;

}


bool SPRSearch::applySPRRound(JointTree &jointTree, int radius, double &bestLoglk, bool blo) {
  return applySPRRoundDigg(jointTree, radius, blo);
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
      blo);
  bestLoglk = jointTree.computeJointLoglk();
  bool foundBetterMove = applyTheBetterMoves(jointTree, 
      betterMoves,
      blo,
      bestLoglk);
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
  jointTree.optimizeParameters(true, false);
  bestLoglk = jointTree.computeJointLoglk();
  while (applySPRRound(jointTree, 5, bestLoglk)) {}
}

