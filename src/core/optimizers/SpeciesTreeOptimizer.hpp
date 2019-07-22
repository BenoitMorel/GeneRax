#pragma once

#include <trees/SpeciesTree.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <string>
#include <maths/DTLRates.hpp>
#include <util/enums.hpp>
#include <families/Families.hpp>
#include <memory>

class SpeciesTreeOptimizer {
public:
  SpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &initialFamilies, 
      RecModel model, 
      const std::string &outputDir,
      const std::string &execPath);

  void rootExhaustiveSearch(bool doOptimizeGeneTrees);
  double sprRound(int radius, bool doOptimizeGeneTrees);
  double hybridSprRound(int radius);
  double sprSearch(int radius, bool doOptimizeGeneTrees);
  void ratesOptimization();
  void advancedRatesOptimization(int radius);
  void optimizeGeneTrees(int radius);
private:
  std::shared_ptr<SpeciesTree> _speciesTree;
  std::shared_ptr<PerCoreGeneTrees> _geneTrees;
  Families _currentFamilies;
  RecModel _model;
  std::string _outputDir;
  std::string _execPath;
  unsigned int _geneTreeIteration;
private:
  double computeReconciliationLikelihood(bool doOptimizeGeneTrees, int geneSPRRadius = 1);
  void rootExhaustiveSearchAux(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model, bool doOptimizeGeneTrees, std::vector<unsigned int> &movesHistory, std::vector<unsigned int> &bestMovesHistory, double &bestLL, unsigned int &visits);
  void saveCurrentSpeciesTree();
};


