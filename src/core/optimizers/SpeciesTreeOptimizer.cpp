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
  _geneTreeIteration(0),
  _perSpeciesRatesOptimization(false)
{
  if (speciesTreeFile == "random") {
    _speciesTree = std::make_unique<SpeciesTree>(initialFamilies);
    _speciesTree->setRates(Parameters(0.1, 0.2, 0.1));
  } else {
    _speciesTree = std::make_unique<SpeciesTree>(speciesTreeFile);
    ratesOptimization();
  }
  _speciesTree->saveToFile(FileSystem::joinPaths(_outputDir, "starting_species_tree.newick"), true);
  std::string subsamplesPath = FileSystem::joinPaths(_outputDir, "subsamples");
  FileSystem::mkdir(FileSystem::joinPaths(_outputDir, "sub_genes_opt"), true);
  FileSystem::mkdir(subsamplesPath, true);
  saveCurrentSpeciesTree();
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
      double ll = computeLikelihood(doOptimizeGeneTrees, 1);
      visits++;
      if (ll > bestLL) {
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
  Logger::timed << "Trying to re-root the species tree" << std::endl;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  double bestLL = computeLikelihood(doOptimizeGeneTrees, 1);
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

double SpeciesTreeOptimizer::sprRound(int radius)
{
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  double bestLL = computeLikelihood(false, 1);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
      double newLL = computeLikelihood(false, 1);
      if (newLL > bestLL) {
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


double SpeciesTreeOptimizer::sortedSprRound(int radius, double bestLL)
{
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  std::vector<EvaluatedMove> evaluatedMoves;
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
      EvaluatedMove em;
      em.prune = prune;
      em.regraft = regraft;
      em.ll = computeLikelihood(false, 1);
      evaluatedMoves.push_back(em);
      SpeciesTreeOperator::reverseSPRMove(*_speciesTree, em.prune, rollback);
    }
  }
  std::sort(evaluatedMoves.begin(), evaluatedMoves.end(), less_than_evaluatedmove());
  unsigned int movesToTry = 50;
  if (movesToTry > evaluatedMoves.size()) {
    movesToTry = evaluatedMoves.size();
  }
  for (unsigned int i = 0; i < movesToTry; ++i) {
    auto &em = evaluatedMoves[i];
    unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, em.prune, em.regraft);
    double newLL = computeLikelihood(true);
    if (newLL > bestLL + 0.001) {
      Logger::info << "Best tree " << _speciesTree->getHash() << " ll=" << newLL << std::endl;
      saveCurrentSpeciesTree();
      optimizeGeneTrees(1, true);
      return newLL;
    }
    SpeciesTreeOperator::reverseSPRMove(*_speciesTree, em.prune, rollback);
  } 
  return bestLL;
}

double SpeciesTreeOptimizer::sprSearch(int radius, bool doOptimizeGeneTrees)
{
  double bestLL = computeLikelihood(doOptimizeGeneTrees, doOptimizeGeneTrees);
  Logger::timed << "Starting species SPR search (" << (doOptimizeGeneTrees ? "SLOW" : "FAST") << ", radius=" << radius << ", bestLL=" << bestLL << ")" <<std::endl;
  double newLL = bestLL;
  //Logger::info << "Initial " << likelihoodName(doOptimizeGeneTrees) << ": " << newLL << std::endl;
  do {
    bestLL = newLL;
    if (doOptimizeGeneTrees) {
      newLL = sortedSprRound(radius, bestLL); 
    } else {
      newLL = sprRound(radius);
    }
  } while (newLL - bestLL > 0.001);
  saveCurrentSpeciesTree();
  return newLL;
}
  
void SpeciesTreeOptimizer::ratesOptimization()
{ 
  if (_perSpeciesRatesOptimization) {
    perSpeciesRatesOptimization();
    return;
  }
  Logger::timed << "Starting DTL rates optimization" << std::endl;
  auto rates = DTLOptimizer::optimizeParametersGlobalDTL(*_geneTrees, _speciesTree->getTree(), _model);
  //Logger::info << " Best rates: " << rates << std::endl;
  _speciesTree->setRates(rates);
}
  
void SpeciesTreeOptimizer::perSpeciesRatesOptimization()
{
  Logger::info << "Starting DTL rates vector optimization" << std::endl;
  Parameters rates = DTLOptimizer::optimizeParametersPerSpecies(*_geneTrees, _speciesTree->getTree(), _model);
  _speciesTree->setRatesVector(rates);
}

void SpeciesTreeOptimizer::saveCurrentSpeciesTree(std::string name)
{
  _speciesTree->saveToFile(FileSystem::joinPaths(_outputDir, name), true);
  ParallelContext::barrier();
}

double SpeciesTreeOptimizer::optimizeGeneTrees(int radius, bool inPlace)
{
  saveCurrentSpeciesTree("proposal_species_tree.newick");
  std::string speciesTree = FileSystem::joinPaths(_outputDir, "proposal_species_tree.newick");
  RecOpt recOpt = Simplex;
  bool perFamilyDTLRates = false;
  bool rootedGeneTree = true;
  double recWeight = 1.0;
  bool useSplitImplem = true;
  long int sumElapsedSPR = 0;
  bool pruneSpeciesTree = false;
  auto rates = _speciesTree->getRates();
  Logger::mute();
  std::string resultName = "proposals";
  Families families = _currentFamilies;
  GeneTreeSearchMaster::optimizeGeneTrees(families, 
      _model, rates, _outputDir, resultName, 
      _execPath, speciesTree, recOpt, perFamilyDTLRates, rootedGeneTree, pruneSpeciesTree, 
      recWeight, true, true, radius, _geneTreeIteration, 
        useSplitImplem, sumElapsedSPR, inPlace);
  _geneTreeIteration++;
  Logger::unmute();
  double totalLibpllLL, totalRecLL;
  Routines::gatherLikelihoods(families, totalLibpllLL, totalRecLL);
  //Logger::info << "optll " << totalLibpllLL + totalRecLL << " " << totalLibpllLL << " " << totalRecLL << std::endl;
  if (inPlace) {
    ratesOptimization();
    Logger::info << "Uptimed the gene trees " << totalLibpllLL + totalRecLL << " " << totalLibpllLL << " " << totalRecLL << std::endl;
  }
  return totalLibpllLL + totalRecLL;
  //_geneTrees = std::make_unique<PerCoreGeneTrees>(_currentFamilies);
}
  
double SpeciesTreeOptimizer::computeLikelihood(bool doOptimizeGeneTrees, int geneSPRRadius)
{
  if (doOptimizeGeneTrees) {
    auto initialFamilies = _currentFamilies;
    double jointLL = optimizeGeneTrees(geneSPRRadius);
    _currentFamilies = initialFamilies;
    _geneTrees = std::make_unique<PerCoreGeneTrees>(_currentFamilies);
    return jointLL;
  }
  return _speciesTree->computeReconciliationLikelihood(*_geneTrees, _model);
}


void SpeciesTreeOptimizer::inferSpeciesTreeFromSamples(unsigned int sampleSize, const std::string &outputSpeciesId)
{
  Logger::info << "Infer species tree from gene tree samples " << outputSpeciesId << std::endl;
  std::string subsamplesPath = FileSystem::joinPaths(_outputDir, "subsamples");
  std::string subOutputDir = FileSystem::joinPaths(subsamplesPath, outputSpeciesId);
  FileSystem::mkdir(subOutputDir, true);
  auto speciesTree = _speciesTree->buildRandomTree();
  std::string speciesTreePath = getSpeciesTreePath(outputSpeciesId);
  ParallelContext::barrier();
  speciesTree->saveToFile(speciesTreePath, true);
  ParallelContext::barrier();
  Families subFamilies;
  for (unsigned int i = 0; i < sampleSize; ++i) {
    subFamilies.push_back(_currentFamilies[rand() % _currentFamilies.size()]);
  }
  SpeciesTreeOptimizer subOptimizer(speciesTreePath, subFamilies, _model, subOutputDir, _execPath);
  for (unsigned int radius = 0; radius < 7; ++radius) {
    subOptimizer.ratesOptimization();
    subOptimizer.sprSearch(radius, false);
    subOptimizer.rootExhaustiveSearch(false);
  }
  subOptimizer.saveCurrentSpeciesTree();
}

std::string SpeciesTreeOptimizer::getSpeciesTreePath(const std::string &speciesId)
{
  std::string subsamplesPath = FileSystem::joinPaths(_outputDir, "subsamples");
  std::string subOutputDir = FileSystem::joinPaths(subsamplesPath, speciesId);
  return FileSystem::joinPaths(subOutputDir, "inferred_species_tree.newick");
}

void SpeciesTreeOptimizer::optimizeGeneTreesFromSamples(const std::unordered_set<std::string> &speciesIds, const std::string &stepId)
{
  Logger::info << "Optimize gene trees from sample species trees " << stepId << std::endl;
  std::string geneOptPath = FileSystem::joinPaths(_outputDir, "sub_genes_opt");
  std::string stepPath = FileSystem::joinPaths(geneOptPath, stepId);
  FileSystem::mkdir(stepPath, true);
  auto families = _currentFamilies;
  std::random_shuffle(families.begin(), families.end());
  unsigned int i = 0;
  unsigned int step = families.size() / speciesIds.size();
  for (auto &speciesId: speciesIds) {
    std::string subStepPath = FileSystem::joinPaths(stepPath, std::string("substep_" + std::to_string(i)));
    FileSystem::mkdir(subStepPath, true);
    unsigned int begin = step * i;
    unsigned int end = std::min<unsigned int>(begin + step, families.size());
    auto subFamilies = Families(families.begin() + begin, families.begin() + end);
    auto speciesTreePath = getSpeciesTreePath(speciesId);
    SpeciesTreeOptimizer subOptimizer(speciesTreePath, subFamilies, _model, subStepPath, _execPath);
    subOptimizer.optimizeGeneTrees(1, true);
    ++i;
  }
  _geneTrees = std::make_unique<PerCoreGeneTrees>(_currentFamilies);
}
