#pragma once

#include <trees/SpeciesTree.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <string>
#include <maths/Parameters.hpp>
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
  void setModel(RecModel model) {_model = model;}
  void setPerSpeciesRatesOptimization(bool perSpecies) {_perSpeciesRatesOptimization = perSpecies;}

  void rootExhaustiveSearch(bool doOptimizeGeneTrees);
  double sprRound(int radius);
  double sortedSprRound(int radius, double bestLL);
  double sprSearch(int radius, bool doOptimizeGeneTrees);
  void ratesOptimization();
  double optimizeGeneTrees(int radius, bool inPlace = false);

  void inferSpeciesTreeFromSamples(unsigned int sampleSize, const std::string &outputSpeciesId);
  void optimizeGeneTreesFromSamples(const std::unordered_set<std::string> &speciesIds, const std::string &stepId);

  double getReconciliationLikelihood() {return computeLikelihood(false);}
  double computeLikelihood(bool doOptimizeGeneTrees, int geneSPRRadius = 1);
  void saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
  
private:
  std::shared_ptr<SpeciesTree> _speciesTree;
  std::shared_ptr<PerCoreGeneTrees> _geneTrees;
  Families _currentFamilies;
  RecModel _model;
  std::string _outputDir;
  std::string _execPath;
  unsigned int _geneTreeIteration;
  bool _perSpeciesRatesOptimization;
private:
  void perSpeciesRatesOptimization();
  void rootExhaustiveSearchAux(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model, bool doOptimizeGeneTrees, std::vector<unsigned int> &movesHistory, std::vector<unsigned int> &bestMovesHistory, double &bestLL, unsigned int &visits);
  std::string getSpeciesTreePath(const std::string &speciesId);
};


