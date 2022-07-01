#include "SearchUtils.hpp"
#include <trees/JointTree.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>

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


void SearchUtils::testMove(JointTree &jointTree,
    SPRMove &move,
    double &newLoglk,
    bool blo,
    std::unordered_map<unsigned int, double> *treeHashScores)
{
  jointTree.applyMove(move); 
  if (treeHashScores) {
    auto hash = jointTree.getUnrootedTreeHash();
    auto it = treeHashScores->find(hash);
    if (it != treeHashScores->end()) {
      newLoglk = it->second;
      jointTree.rollbackLastMove();
      move.setScore(newLoglk);
      return;
    }
  }
  double recLoglk = jointTree.computeReconciliationLoglk();
  if (blo) {
    jointTree.optimizeMove(move);
  }
  newLoglk = recLoglk +  jointTree.computeLibpllLoglk(false);
  move.setScore(newLoglk);
  jointTree.rollbackLastMove();
  if (treeHashScores) {
    treeHashScores->insert({jointTree.getUnrootedTreeHash(), newLoglk});
  } 
}


struct less_than_move_ptr
{
  inline bool operator() (const std::shared_ptr<SPRMove> &m1, 
      const std::shared_ptr<SPRMove> &m2)
  {
    if (m1->getScore() != m2->getScore()) {
      return (m1->getScore() < m2->getScore());
    } else if (m1->getPruneIndex() != m2->getPruneIndex()) {
      return m1->getPruneIndex() <  m2->getPruneIndex();
    } else {
      return m1->getRegraftIndex() < m2->getRegraftIndex();
    }
  }
};


bool SearchUtils::findBetterMoves(JointTree &jointTree,
  std::vector<std::shared_ptr<SPRMove> > &allMoves,
  std::vector<std::shared_ptr<SPRMove> > &sortedBetterMoves,
  bool blo,
  bool check)
{
  double initialReconciliationLoglk = jointTree.computeReconciliationLoglk();
  double initialLibpllLoglk = jointTree.computeLibpllLoglk();
  double initialLoglk = initialLibpllLoglk + initialReconciliationLoglk;
  auto begin = ParallelContext::getBegin(static_cast<unsigned int>(allMoves.size()));
  auto end = ParallelContext::getEnd(static_cast<unsigned int>(allMoves.size()));
  std::vector<unsigned int> localGoodIndices;
  std::vector<double> localGoodLL;
  std::vector<unsigned int> goodIndices;
  std::vector<double> goodLL;
  double epsilon = fabs(initialLoglk) > 10000.0 ? 0.5 : 0.001; 
  for (auto i = begin; i < end; ++i) {
    auto loglk = initialLoglk;
    SearchUtils::testMove(jointTree, *allMoves[i], 
        loglk,
        blo);
    if (loglk - initialLoglk > epsilon) {
      localGoodIndices.push_back(i);
      localGoodLL.push_back(allMoves[i]->getScore());
      assert(allMoves[i]->getScore() != 0.0);
    }
  }
  ParallelContext::concatenateHetherogeneousUIntVectors(localGoodIndices, 
      goodIndices); 
  ParallelContext::concatenateHetherogeneousDoubleVectors(localGoodLL, 
      goodLL);
  assert(ParallelContext::isIntEqual(goodIndices.size()));
  assert(goodIndices.size() == goodLL.size());
  for (unsigned int i = 0; i < goodIndices.size(); ++i) {
    unsigned int index = goodIndices[i];
    double ll = goodLL[i];
    assert(ParallelContext::isIntEqual(index));
    assert(ParallelContext::isDoubleEqual(ll));
    auto move = allMoves[index];
    move->setScore(ll);
    sortedBetterMoves.push_back(move);
  }
  for (auto move: sortedBetterMoves) {
    assert(ParallelContext::isIntEqual(move->getPruneIndex()));
  }
  std::sort(sortedBetterMoves.rbegin(), sortedBetterMoves.rend(), less_than_move_ptr()); 
  for (auto move: sortedBetterMoves) {
    assert(ParallelContext::isIntEqual(move->getPruneIndex()));
  }
  jointTree.getReconciliationEvaluation().invalidateAllCLVs();
  return 0 < sortedBetterMoves.size();
}


static void diggRecursive(JointTree &jointTree,
    std::unordered_map<unsigned int, double> &treeHashScores,
    DiggStopper &stopper,
    corax_unode_t *pruneNode,
    corax_unode_t *regraftNode,
    std::vector<unsigned int> &path,
    unsigned int radius,
    unsigned int maxRadius,
    unsigned int additionalRadius,
    bool blo,
    double &bestLL,
    double &bestLLAmongPrune,
    SPRMove &bestMove)
{
  // should we stop?
  if (radius >= maxRadius) {
    return;
  }
  // test the current move
  SPRMove move(pruneNode->node_index,
      regraftNode->node_index,
      path);
  double diff = 1.0;
  if (isValidSPRMove(pruneNode, regraftNode)) {
    double newLL = 0.0;
    SearchUtils::testMove(jointTree, 
      move,
      newLL,
      blo);
    diff = newLL - bestLL;
    if (newLL > bestLLAmongPrune) {
      bestMove = move;
      bestLLAmongPrune = newLL;
      if (newLL > bestLL) {
        Logger::info << "better " << newLL << " prune=" << pruneNode->node_index << " r=" << radius <<  std::endl;
        bestLL = newLL;
        maxRadius += additionalRadius;
      }
    }
    stopper.treatDiff(diff, radius);
  }

  // continue 
  radius += 1;
  if (regraftNode->next && radius < maxRadius && !stopper.doStop(diff, radius)) {
    auto left = regraftNode->next->back;
    auto right = regraftNode->next->next->back;
    path.push_back(regraftNode->node_index);
    diggRecursive(jointTree, treeHashScores, stopper, pruneNode, left, path, radius, maxRadius,
        additionalRadius, blo, bestLL, bestLLAmongPrune, bestMove);
    diggRecursive(jointTree, treeHashScores, stopper, pruneNode, right, path, radius, maxRadius,
        additionalRadius, blo, bestLL, bestLLAmongPrune, bestMove);
    path.pop_back();
  }
}

bool SearchUtils::diggBestMoveFromPrune(JointTree &jointTree,
    std::unordered_map<unsigned int, double> &treeHashScores,
    DiggStopper &stopper,
    unsigned int pruneIndex,
    unsigned int maxRadius,
    unsigned int additionalRadius,
    double &bestLL,
    double &bestLLAmongPrune,
    bool blo,
    SPRMove &bestMove)
{
  auto pruneNode = jointTree.getNode(pruneIndex);
  assert(pruneNode->next);
  auto regraft1 = pruneNode->next->back;
  auto regraft2 = pruneNode->next->next->back;
  std::vector<unsigned int> path;
  bestLLAmongPrune = -99999999999.0;
  double initialLL = bestLL;
  diggRecursive(jointTree, treeHashScores, stopper, pruneNode, regraft1, path, 
      0, maxRadius, additionalRadius, blo, 
      bestLL, bestLLAmongPrune, bestMove);
  diggRecursive(jointTree, treeHashScores, stopper, pruneNode, regraft2, path, 
      0, maxRadius, additionalRadius, blo, 
      bestLL, bestLLAmongPrune, bestMove);
  return  (bestLLAmongPrune - initialLL) > 0.1;
}



