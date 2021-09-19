#include "UNNISearch.hpp"  

#include <IO/Logger.hpp>

static void uswap(pll_unode_t *A,
    pll_unode_t *B) 
{
  auto temp = A->back;
  A->back = B->back;
  B->back = temp;
  A->back->back = A;
  B->back->back = B;
}

void UNNIMove::apply()
{
  uswap(getB(), getC());
}


double UNNISearch::runRound()
{
  double bestScore = _evaluator.eval(_tree);
  Logger::timed << "Starting NNI round, score=" << bestScore << std::endl;
  UNNIMove bestMove(nullptr, false);
  std::vector<bool> bools;
  bools.push_back(true);
  bools.push_back(false);
  for (auto edge: _tree.getBranchesDeterministic()) {
    if (!edge->next or !edge->back->next) {
      continue;
    }
    for (bool left: bools) {
      UNNIMove move(edge, left);
      double newScore = _evaluator.evalNNI(_tree, 
          move);
      if (newScore > bestScore) {
        bestScore = newScore;
        bestMove = move;
        Logger::info << "Better tree, score=" << bestScore << std::endl;
      }
    }
  }
  if (bestMove.edge != nullptr) {
    bestMove.apply();
  }
  return bestScore;
}

double UNNISearch::search()
{
  double bestScore = _evaluator.eval(_tree);
  double epsilon = 0.0;
  do {
    double newScore = runRound();
    epsilon = newScore - bestScore;
    bestScore = newScore;
  } while (epsilon > 0.000000001);
  return bestScore;
}


