#include "DatedSpeciesTreeSearch.hpp"
#include <IO/Logger.hpp>
#include <trees/SpeciesTree.hpp>
#include <optimizers/SpeciesTreeOptimizer.hpp>
#include <search/SpeciesTransferSearch.hpp>
  

static void saveDatedTree(SpeciesTree &speciesTree,
  const SpeciesSearchState &state)
{
  if (state.pathToBestSpeciesTree.size() == 0) {
    return;
  }
  speciesTree.getDatedTree().rescaleBranchLengths();
  speciesTree.saveToFile(state.pathToBestSpeciesTree, true);
  ParallelContext::barrier();
}


static double optimizeDatesNaive(SpeciesTree &speciesTree,
    SpeciesSearchState &state,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation)
{
  auto &tree = speciesTree.getDatedTree();
  unsigned int max = tree.getOrderedSpeciations().size();
  auto bestLL = evaluation.computeLikelihood();
  bool tryAgain = true;
  Logger::timed << "Starting a new naive dating search from ll=" << bestLL << std::endl;
  while (tryAgain) {
    tryAgain = false;
    for (unsigned int rank = 0; rank < max; ++rank) {
      if (!tree.moveUp(rank)) {
        continue;
      }
      evaluation.onSpeciesDatesChange();
      auto ll = evaluation.computeLikelihood();
      if (ll > state.bestLL) {
        // best tree over all search iterations, we save it!
        state.bestLL = ll;
        saveDatedTree(speciesTree, state);
      }
      if (ll > bestLL) {
        // best tree for this search iteration
        // (we might have found a better likelihood
        // in another iteration)
        tryAgain = true;
        bestLL = ll;
        rank -= std::min((unsigned int)2, rank);
      } else {
        tree.moveUp(rank);
      }
      assert(tree.isConsistent());
    }
    Logger::timed << " end of round, ll= " << bestLL <<  std::endl;
  }
  evaluation.onSpeciesDatesChange();
  Logger::timed << "End of naive dating search, ll= " << bestLL <<  std::endl;
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
  auto maxDisplacement = size_t(sqrt(double(N)) * 2.0 * perturbation);
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
      SpeciesSearchState &state,
      double currentLL,
      bool thorough)
{
  auto bestLL = currentLL;
  if (!evaluation.isDated()) {
    return bestLL;
  }
  Logger::timed << "Optimizing dates, ll=" << bestLL << std::endl;
  optimizeDatesNaive(speciesTree, state, evaluation); 
  bestLL = evaluation.computeLikelihood();
  unsigned int unsuccessfulTrials = 0;
  const unsigned int maxTrials = 2;
  auto &tree = speciesTree.getDatedTree();
  while (thorough && unsuccessfulTrials < maxTrials) {
    auto backup = tree.getBackup();
    double perturbation = double(unsuccessfulTrials + 1) / double(maxTrials);
    perturbateDates(speciesTree, perturbation);
    evaluation.onSpeciesDatesChange();
    auto ll = evaluation.computeLikelihood();
    optimizeDatesNaive(speciesTree, state, evaluation);
    ll = evaluation.computeLikelihood();
    if (ll <= bestLL) {
      tree.restore(backup);
      unsuccessfulTrials++;
    } else {
      unsuccessfulTrials = 0;
      bestLL = ll;
      Logger::timed << " better ll=" << bestLL << std::endl;
    }
  } 
  Logger::timed << "after optimizating dates, ll = " << bestLL << std::endl;
  return bestLL;
}


unsigned int getTransferScore(SpeciesTree &speciesTree,
  const TransferFrequencies &frequencies)
{
  std::unordered_map<std::string, unsigned int> labelsToIds;
  speciesTree.getLabelsToId(labelsToIds);
  unsigned int score = 0;
  for (unsigned int from = 0; from < frequencies.count.size(); ++from) {
    for (unsigned int to = 0; to < frequencies.count.size(); ++to) {
      auto count = frequencies.count[from][to];
      if (!count) {
        continue;
      }
      auto src = labelsToIds[frequencies.idToLabel[from]];
      auto dest = labelsToIds[frequencies.idToLabel[to]];
      if (speciesTree.getDatedTree().canTransferUnderRelDated(src, dest)) {
        score += count;
      }
    }
  }
  return score;
}

class FakeEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
 private:
  SpeciesTree &_speciesTree;
  TransferFrequencies &_frequencies;
public:
  FakeEvaluator(SpeciesTree &speciesTree, 
      TransferFrequencies &frequencies):
    _speciesTree(speciesTree),
    _frequencies(frequencies)
  {}
  
  virtual double computeLikelihood() {
    return computeLikelihoodFast();
  }
  virtual double computeLikelihoodFast() {
    return getTransferScore(_speciesTree, _frequencies);
  }
  virtual bool providesFastLikelihoodImpl() const {return false;}
  virtual bool isDated() const {return true;}
  virtual double optimizeModelRates(bool) {
    return computeLikelihood();
  }
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void getTransferInformation(SpeciesTree &,
    TransferFrequencies &,
    PerSpeciesEvents &,
    PerCorePotentialTransfers &) {
    assert(false); 
  }
  virtual void fillPerFamilyLikelihoods(PerFamLL &) {
    assert(false);
  }
  virtual bool pruneSpeciesTree() const {return false;}
};

double DatedSpeciesTreeSearch::optimizeDatesFromReconciliation(SpeciesTree &speciesTree,
    SpeciesSearchState &state,  
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation)
{
  Logger::timed << "Optimize dates from reconciliations, bestLL=" << evaluation.computeLikelihood() << std::endl;
  TransferFrequencies frequencies;
  PerSpeciesEvents perSpeciesEvents;
  PerCorePotentialTransfers potentialTransfers;
  evaluation.getTransferInformation(speciesTree, frequencies,
      perSpeciesEvents, potentialTransfers);
  FakeEvaluator fakeEvaluator(speciesTree, frequencies); 
  perturbateDates(speciesTree, 100.0);
  optimizeDatesNaive(speciesTree,
      state,
      fakeEvaluator);      
  
  auto bestScore = fakeEvaluator.computeLikelihood();
  unsigned int unsuccessfulTrials = 0;
  const unsigned int maxTrials = 10;
  auto &tree = speciesTree.getDatedTree();
  while (unsuccessfulTrials < maxTrials) {
    auto backup = tree.getBackup();
    double perturbation = double(unsuccessfulTrials + 1) / double(maxTrials);
    perturbateDates(speciesTree, perturbation);
    optimizeDatesNaive(speciesTree, state, fakeEvaluator);
    auto score = fakeEvaluator.computeLikelihood();
    score = fakeEvaluator.computeLikelihood();
    if (score <= bestScore) {
      tree.restore(backup);
      unsuccessfulTrials++;
    } else {
      unsuccessfulTrials = 0;
      bestScore = score;
      Logger::timed << " better score=" << bestScore  << std::endl;//" ll=" << evaluation.computeLikelihood()<< std::endl;
    }
  }
  evaluation.onSpeciesDatesChange();
  
  return evaluation.computeLikelihood();
}
