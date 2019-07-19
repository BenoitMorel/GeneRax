#include "SpeciesTreeOptimizer.hpp"

#include <optimizers/DTLOptimizer.hpp>
#include <IO/FileSystem.hpp>

SpeciesTreeOptimizer::SpeciesTreeOptimizer(const std::string speciesTreeFile, 
    const Families &initialFamilies, 
    RecModel model, 
    const std::string &outputDir):
  _speciesTree(0),
  _geneTrees(std::make_unique<PerCoreGeneTrees>(initialFamilies)),
  _currentFamilies(initialFamilies),
  _model(model),
  _outputDir(outputDir)
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




  
void rootExhaustiveSearchAux(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model, std::vector<unsigned int> &movesHistory, std::vector<unsigned int> &bestMovesHistory, double &bestLL, unsigned int &visits)
{
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  for (auto direction: moves) {
    if (SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      movesHistory.push_back(direction);
      SpeciesTreeOperator::changeRoot(speciesTree, direction);
      double ll = speciesTree.computeReconciliationLikelihood(geneTrees, model);
      visits++;
      if (ll > bestLL) {
        bestLL = ll;
        bestMovesHistory = movesHistory;
      }
      rootExhaustiveSearchAux(speciesTree, geneTrees, model, movesHistory, bestMovesHistory, bestLL, visits);
      SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
      movesHistory.pop_back();
    }
  }
}

void SpeciesTreeOptimizer::rootExhaustiveSearch()
{
  Logger::info << "Trying to re-root the species tree" << std::endl;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  double bestLL = _speciesTree->computeReconciliationLikelihood(*_geneTrees, _model);
  unsigned int visits = 1;
  movesHistory.push_back(0);
  rootExhaustiveSearchAux(*_speciesTree, *_geneTrees, _model, movesHistory, bestMovesHistory, bestLL, visits); 
  movesHistory[0] = 1;
  rootExhaustiveSearchAux(*_speciesTree, *_geneTrees, _model, movesHistory, bestMovesHistory, bestLL, visits); 
  assert (visits == 2 * _speciesTree->getTaxaNumber() - 3);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(*_speciesTree, bestMovesHistory[i]);
  }
}
  
double SpeciesTreeOptimizer::sprRound(int radius)
{
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  double bestLL = _speciesTree->computeReconciliationLikelihood(*_geneTrees, _model);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
      double newLL = _speciesTree->computeReconciliationLikelihood(*_geneTrees, _model);
      if (newLL > bestLL) {
        return newLL;
      }
      SpeciesTreeOperator::reverseSPRMove(*_speciesTree, prune, rollback);
    }
  }
  return bestLL;
}

double SpeciesTreeOptimizer::sprSearch(int radius)
{
  Logger::info << "Starting species SPR round (radius=" << radius << ")" <<std::endl;
  double bestLL = _speciesTree->computeReconciliationLikelihood(*_geneTrees, _model);
  double newLL = bestLL;
  do {
    bestLL = newLL;
    Logger::info << "LL = " << bestLL << std::endl;
    newLL = sprRound(radius);
  } while (newLL - bestLL > 0.001);
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
