#include "SpeciesTreeOptimizer.hpp"

#include <optimizers/DTLOptimizer.hpp>
#include <IO/FileSystem.hpp>
#include <routines/GeneTreeSearchMaster.hpp>
#include <routines/Routines.hpp>
#include <algorithm>

SpeciesTreeOptimizer::SpeciesTreeOptimizer(const std::string speciesTreeFile, 
    const Families &initialFamilies, 
    RecModel model, 
    const std::string &outputDir,
    const std::string &execPath):
  _speciesTree(0),
  _geneTrees(std::make_unique<PerCoreGeneTrees>(initialFamilies)),
  _currentFamilies(initialFamilies),
  _model(model),
  _outputDir(outputDir),
  _execPath(execPath),
  _geneTreeIteration(0)
{
  if (speciesTreeFile == "random") {
    _speciesTree = std::make_unique<SpeciesTree>(initialFamilies);
    _speciesTree->setRates(DTLRates(0.1, 0.2, 0.1));
  } else {
    _speciesTree = std::make_unique<SpeciesTree>(speciesTreeFile);
    ratesOptimization();
  }
  _speciesTree->saveToFile(FileSystem::joinPaths(_outputDir, "starting_species_tree.newick"));
}




  
void SpeciesTreeOptimizer::rootExhaustiveSearchAux(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model, bool doOptimizeGeneTrees, std::vector<unsigned int> &movesHistory, std::vector<unsigned int> &bestMovesHistory, double &bestLL, unsigned int &visits)
{
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  for (auto direction: moves) {
    if (SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      movesHistory.push_back(direction);
      SpeciesTreeOperator::changeRoot(speciesTree, direction);
      double ll = computeReconciliationLikelihood(doOptimizeGeneTrees);
      visits++;
      if (ll > bestLL) {
        Logger::info << "Found best root " << std::endl;
        bestLL = ll;
        bestMovesHistory = movesHistory;
      }
      rootExhaustiveSearchAux(speciesTree, geneTrees, model, doOptimizeGeneTrees, movesHistory, bestMovesHistory, bestLL, visits);
      SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
      movesHistory.pop_back();
    }
  }
}

void SpeciesTreeOptimizer::rootExhaustiveSearch(bool doOptimizeGeneTrees)
{
  Logger::info << "Trying to re-root the species tree" << std::endl;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  double bestLL = computeReconciliationLikelihood(doOptimizeGeneTrees);
  unsigned int visits = 1;
  movesHistory.push_back(0);
  rootExhaustiveSearchAux(*_speciesTree, *_geneTrees, _model, doOptimizeGeneTrees, movesHistory, bestMovesHistory, bestLL, visits); 
  movesHistory[0] = 1;
  rootExhaustiveSearchAux(*_speciesTree, *_geneTrees, _model, doOptimizeGeneTrees, movesHistory, bestMovesHistory, bestLL, visits); 
  assert (visits == 2 * _speciesTree->getTaxaNumber() - 3);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(*_speciesTree, bestMovesHistory[i]);
  }
}
  
std::string likelihoodName(bool doOptimizeGeneTrees) {
  return doOptimizeGeneTrees ? std::string("Joint LL") : std::string("Rec LL");
}

double SpeciesTreeOptimizer::sprRound(int radius, bool doOptimizeGeneTrees)
{
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  double bestLL = computeReconciliationLikelihood(doOptimizeGeneTrees, 1);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
      double newLL = computeReconciliationLikelihood(doOptimizeGeneTrees);
      if (newLL > bestLL) {
        Logger::info << "New best " << likelihoodName(doOptimizeGeneTrees) << ": " << newLL << std::endl; 
        return newLL;
      }
      SpeciesTreeOperator::reverseSPRMove(*_speciesTree, prune, rollback);
    }
  }
  return bestLL;
}

struct EvaluatedMove {
  unsigned int prune;
  unsigned int regraft;
  double ll;
};

struct less_than_evaluatedmove
{
  inline bool operator() (const EvaluatedMove& e1, const EvaluatedMove& e2)
  {
    return e1.ll > e2.ll;
  }
};


double SpeciesTreeOptimizer::hybridSprRound(int radius)
{
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  double bestLL = computeReconciliationLikelihood(false, 1);
  std::vector<EvaluatedMove> evaluatedMoves;
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
      EvaluatedMove em;
      em.prune = prune;
      em.regraft = regraft;
      em.ll = computeReconciliationLikelihood(false);
      evaluatedMoves.push_back(em);
      SpeciesTreeOperator::reverseSPRMove(*_speciesTree, em.prune, rollback);
    }
  }
  std::sort(evaluatedMoves.begin(), evaluatedMoves.end(), less_than_evaluatedmove());
  unsigned int movesToTry = 20;
  if (movesToTry > evaluatedMoves.size()) {
    movesToTry = evaluatedMoves.size();
  }
  bestLL = computeReconciliationLikelihood(true, 1);
  srand(42);
  Logger::info << "Mistgabel!!! " << bestLL << std::endl;
  srand(42);
  Logger::info << "Mistgabel!!! " << computeReconciliationLikelihood(true, 1) << std::endl;
  srand(42);
  Logger::info << "Mistgabel!!! " << computeReconciliationLikelihood(true, 1) << std::endl;
  srand(42);
  Logger::info << "Mistgabel!!! " << computeReconciliationLikelihood(true, 1) << std::endl;
  for (unsigned int i = 0; i < movesToTry; ++i) {
    auto &em = evaluatedMoves[i];
    //Logger::info << "evaluated move with initial rec ll = " << em.ll << std::endl;
    unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, em.prune, em.regraft);
    double newLL = computeReconciliationLikelihood(true);
    if (newLL > bestLL + 0.001) {
      Logger::info << "New best " << likelihoodName(true) << ": " << newLL << std::endl;
      optimizeGeneTrees(1);
      return newLL;
    }
    SpeciesTreeOperator::reverseSPRMove(*_speciesTree, em.prune, rollback);
  } 
  return bestLL;
}

double SpeciesTreeOptimizer::sprSearch(int radius, bool doOptimizeGeneTrees)
{
  Logger::info << "Starting species SPR search (" << (doOptimizeGeneTrees ? "SLOW" : "FAST") << ", radius=" << radius << ")" <<std::endl;
  double bestLL = computeReconciliationLikelihood(doOptimizeGeneTrees);
  double newLL = bestLL;
  Logger::info << "Initial " << likelihoodName(doOptimizeGeneTrees) << ": " << newLL << std::endl;
  do {
    bestLL = newLL;
    if (doOptimizeGeneTrees) {
      newLL = hybridSprRound(radius); 
    } else {
      newLL = sprRound(radius, doOptimizeGeneTrees);
    }
  } while (newLL - bestLL > 0.001);
  saveCurrentSpeciesTree();
  Logger::info <<"End of SPR Search" << std::endl;
  return newLL;
}
  
void SpeciesTreeOptimizer::ratesOptimization()
{
  Logger::info << "Starting DTL rates optimization" << std::endl;
  DTLRates rates = DTLOptimizer::optimizeDTLRates(*_geneTrees, _speciesTree->getTree(), _model);
  Logger::info << " Best rates: " << rates << std::endl;
  _speciesTree->setRates(rates);
}

void SpeciesTreeOptimizer::saveCurrentSpeciesTree()
{
  _speciesTree->saveToFile(FileSystem::joinPaths(_outputDir, "inferred_species_tree.newick"));
}

void SpeciesTreeOptimizer::optimizeGeneTrees(int radius)
{
  saveCurrentSpeciesTree();
  std::string speciesTree = FileSystem::joinPaths(_outputDir, "inferred_species_tree.newick");
  RecOpt recOpt = Simplex;
  bool rootedGeneTree = true;
  double recWeight = 1.0;
  bool useSplitImplem = true;
  long int sumElapsedSPR = 0;
  auto rates = _speciesTree->getRates();
  DTLRatesVector ratesVector(rates);
  Logger::mute();
  std::string resultName = "proposals";
  GeneTreeSearchMaster::optimizeGeneTrees(_currentFamilies, 
      _model, ratesVector, _outputDir, resultName, 
      _execPath, speciesTree, recOpt, rootedGeneTree, 
      recWeight, true, radius, _geneTreeIteration, 
        useSplitImplem, sumElapsedSPR);
  _geneTreeIteration++;
  Logger::unmute();
  _geneTrees = std::make_unique<PerCoreGeneTrees>(_currentFamilies);
}
  
void SpeciesTreeOptimizer::advancedRatesOptimization(int radius)
{
  auto initialFamilies = _currentFamilies;
  optimizeGeneTrees(radius);
  ratesOptimization();
  _currentFamilies = initialFamilies;
  _geneTrees = std::make_unique<PerCoreGeneTrees>(_currentFamilies);

}
  
double SpeciesTreeOptimizer::computeReconciliationLikelihood(bool doOptimizeGeneTrees, int geneSPRRadius)
{
  if (doOptimizeGeneTrees) {
    auto initialFamilies = _currentFamilies;
    optimizeGeneTrees(geneSPRRadius);
    double totalLibpllLL = 0.0;
    double totalRecLL = 0.0;
    Logger::mute();
    Routines::gatherLikelihoods(_currentFamilies, totalLibpllLL, totalRecLL);
    Logger::unmute();
    _currentFamilies = initialFamilies;
    _geneTrees = std::make_unique<PerCoreGeneTrees>(_currentFamilies);
    //Logger::info << optll " << totalLibpllLL + totalRecLL << " " << totalLibpllLL << " " << totalRecLL << std::endl;
    return totalLibpllLL + totalRecLL;
  }
  return _speciesTree->computeReconciliationLikelihood(*_geneTrees, _model);
}
