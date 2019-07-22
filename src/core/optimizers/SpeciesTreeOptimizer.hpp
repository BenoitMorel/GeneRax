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
      const std::string &outputDir);

  void rootExhaustiveSearch();
  double sprRound(int radius);
  double sprSearch(int radius);
  void ratesOptimization();
  void optimizeGeneTrees(int radius, const std::string &execPath);
private:
  std::shared_ptr<SpeciesTree> _speciesTree;
  std::shared_ptr<PerCoreGeneTrees> _geneTrees;
  Families _currentFamilies;
  RecModel _model;
  std::string _outputDir;

private:
  void saveCurrentSpeciesTree();
};


