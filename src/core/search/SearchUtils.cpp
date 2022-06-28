#include "SearchUtils.hpp"
#include <trees/JointTree.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>



void SearchUtils::testMove(JointTree &jointTree,
    SPRMove &move,
    double initialReconciliationLoglk,
    double initialLibpllLoglk,
#ifdef SUPER_OPTIM
    double &averageReconciliationDiff,
#else
    double &,
#endif
    double &newLoglk,
    bool blo,
    bool check
    )
{
  double initialLoglk = initialReconciliationLoglk + initialLibpllLoglk;
  jointTree.applyMove(move); 
  double recLoglk = jointTree.computeReconciliationLoglk();
#ifdef SUPER_OPTIM
  double improvement = recLoglk - initialReconciliationLoglk;
  averageReconciliationDiff *= 50;
  averageReconciliationDiff += improvement;
  averageReconciliationDiff /= 51;
  if (improvement < averageReconciliationDiff) {
    jointTree.rollbackLastMove();
    if(check && fabs(initialLoglk - jointTree.computeJointLoglk()) > 0.000001) {
      std::cerr.precision(17);
      std::cerr << "small rollback lead to different likelihoods: " << initialLoglk
        << " " << jointTree.computeJointLoglk() << std::endl;
      std::cerr << " rank " << ParallelContext::getRank() << std::endl;
      exit(1);
    }
    return;
  }
#endif
  if (blo) {
    jointTree.optimizeMove(move);
  }
  newLoglk = recLoglk +  jointTree.computeLibpllLoglk(false);
  move.setScore(newLoglk);
  jointTree.rollbackLastMove();
  if(check) {
    auto rbLoglk = jointTree.computeJointLoglk();
    if (fabs(initialLoglk - rbLoglk) > 0.000001) {
      jointTree.printLoglk();
      std::cerr.precision(17);
      std::cerr << "rollback lead to different likelihoods: " << initialLoglk
        << " " << rbLoglk << std::endl;
      std::cerr << "recomputing the ll again: " << jointTree.computeJointLoglk() << std::endl;
      std::cerr << " rank " << ParallelContext::getRank() << std::endl;
      exit(1);
    }
  }
}

//#define STOP
bool SearchUtils::findBestMove(JointTree &jointTree,
    std::vector<std::shared_ptr<SPRMove> > &allMoves,
    double &bestLoglk,
    unsigned int &bestMoveIndex,
    bool blo,
    bool check)
{
  bestMoveIndex = static_cast<unsigned int>(-1);
  double initialLoglk = bestLoglk; //jointTree.computeJointLoglk();
  double initialReconciliationLoglk = jointTree.computeReconciliationLoglk();
  double initialLibpllLoglk = jointTree.computeLibpllLoglk();
  double averageReconciliationDiff = 0;
  double error = fabs(initialLoglk - 
      (initialReconciliationLoglk + initialLibpllLoglk));
  if (error > 0.01)
  {
    Logger::info << "WARNING: potential numerical issue in SearchUtils::findBestMove " << error << std::endl;
  }
  auto begin = ParallelContext::getBegin(static_cast<unsigned int>(allMoves.size()));
  auto end = ParallelContext::getEnd(static_cast<unsigned int>(allMoves.size()));
  unsigned int bestRank = 0;
#ifdef STOP
  // ensure that all cores recieve the same number of tasks, 
  // to avoid deadlock when synchronizing
  while ((end - begin) * ParallelContext::getSize() < allMoves.size()) {
    if (begin > 0) {
      begin -= 1;
    } else {
      assert(end < allMoves.size());
      end += 1;
    }
  }
#endif
  for (auto i = begin; i < end; ++i) {
    auto loglk = bestLoglk;
    SearchUtils::testMove(jointTree, *allMoves[i], 
        initialReconciliationLoglk,
        initialLibpllLoglk, 
        averageReconciliationDiff,
        loglk,
        blo,
        check);
    if (loglk > bestLoglk) {
      bestLoglk = loglk;
      bestMoveIndex = i;
    }
#ifdef STOP
    if ((begin - i) % 1 == 0) {
      ParallelContext::getMax(bestLoglk, bestRank);
      ParallelContext::broadcastUInt(bestRank, bestMoveIndex);
      if (bestMoveIndex != static_cast<unsigned int>(-1)) {
        return true;
      }
    }
#endif
  }
  ParallelContext::getMax(bestLoglk, bestRank);
  ParallelContext::broadcastUInt(bestRank, bestMoveIndex);
  Logger::info << "best;; " << bestLoglk << " " << bestRank << std::endl;
  jointTree.getReconciliationEvaluation().invalidateAllCLVs();
  return bestMoveIndex != static_cast<unsigned int>(-1);
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
  double averageReconciliationDiff = 0;
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
        initialReconciliationLoglk,
        initialLibpllLoglk, 
        averageReconciliationDiff,
        loglk,
        blo,
        check);
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


