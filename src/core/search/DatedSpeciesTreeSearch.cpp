#include "DatedSpeciesTreeSearch.hpp"
#include <IO/Logger.hpp>
#include <trees/SpeciesTree.hpp>
#include <optimizers/SpeciesTreeOptimizer.hpp>

static double optimizeDatesNaive(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation)
{
  auto &tree = speciesTree.getDatedTree();
  unsigned int max = tree.getOrderedSpeciations().size();
  auto bestLL = evaluation.computeLikelihood();
  bool tryAgain = true;
  while (tryAgain) {
    tryAgain = false;
    //Logger::timed << "new naive dated round ll= " << bestLL <<  std::endl;
    for (unsigned int rank = 0; rank < max; ++rank) {
      if (!tree.moveUp(rank)) {
        continue;
      }
      auto ll = evaluation.computeLikelihood();
      if (ll > bestLL) {
        tryAgain = true;
        bestLL = ll;
        rank -= std::min((unsigned int)2, rank);
      } else {
        tree.moveUp(rank);
      }
      assert(tree.isConsistent());
    }
  }
  //Logger::timed << "end naive dated round ll= " << bestLL <<  std::endl;
  return bestLL;
}


// Randomly perturbate the order of speciation events in the dated 
// species tree. The level of perturbation is proportional with 
// perturbation, which is typicall between 0 and 1 (but can be greater)
static void perturbateDates(SpeciesTree &speciesTree, double perturbation = 1.0)
{
  auto &tree = speciesTree.getDatedTree();
  size_t N = tree.getOrderedSpeciations().size();
  auto perturbations = size_t(double((2 * N) * perturbation));
  auto maxDisplacement = size_t(double(N / 3) * perturbation);
  if (maxDisplacement < 2) {
    maxDisplacement = 2;
  }
  for (size_t i = 0; i < perturbations; ++i) {
    auto rank = Random::getInt() % N;
    auto  displacement =  1 + (Random::getInt() % maxDisplacement);
    auto nodesToMove = static_cast<size_t>((Random::getInt() % 2) + 1);
    for (size_t k = 0; k < nodesToMove; ++k) {
      for (size_t j = 0; j < displacement; ++j) {
        tree.moveUp(rank + k - j);
      }
    }
  }
}

double DatedSpeciesTreeSearch::optimizeDates(
      SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
      SpeciesSearchState &,
      double currentLL,
      bool thorough)
{
  auto bestLL = currentLL;
  if (!evaluation.isDated()) {
    return bestLL;
  }
  Logger::timed << "Optimizing dates, ll=" << bestLL << std::endl;
  optimizeDatesNaive(speciesTree, evaluation); 
  bestLL = evaluation.computeLikelihood();
 
  unsigned int unsuccessfulTrials = 0;
  const unsigned int maxTrials = 3;
  auto &tree = speciesTree.getDatedTree();
  while (thorough && unsuccessfulTrials < maxTrials) {
    auto backup = tree.getBackup();
    double perturbation = double(unsuccessfulTrials + 1) / double(maxTrials);
    perturbateDates(speciesTree, perturbation);
    auto ll = evaluation.computeLikelihood();
    optimizeDatesNaive(speciesTree, evaluation);
    ll = evaluation.computeLikelihood();
    if (ll <= bestLL) {
      tree.restore(backup);
      unsuccessfulTrials++;
    } else {
      unsuccessfulTrials = 0;
      bestLL = ll;
    }
  } 
  Logger::timed << "after optimizating dates, ll = " << bestLL << std::endl;
  return bestLL;
}

  
