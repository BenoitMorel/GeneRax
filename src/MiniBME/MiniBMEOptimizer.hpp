#pragma once

#include <search/SpeciesRootSearch.hpp>
#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>
#include <vector>
#include <maths/ModelParameters.hpp>
#include <NJ/MiniBME.hpp>


class MiniBMEEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  MiniBMEEvaluator(PLLRootedTree &rootedTree,
      const Families &families,
      bool missingData);
  virtual ~MiniBMEEvaluator() {}
  virtual double computeLikelihood();
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const {return false;}
  virtual double optimizeModelRates(bool thorough = false) {
    return computeLikelihood();
  }
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL);
  virtual void getTransferInformation(PLLRootedTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents);
  virtual bool pruneSpeciesTree() const {return false;}
private:
  PLLRootedTree *_speciesTree;
  const Families *_families;
  MiniBME _miniBME;
};


class MiniBMEOptimizer: public SpeciesTree::Listener {
public:
  MiniBMEOptimizer(const std::string speciesTreeFile, 
      const Families &families,
      bool missingData,
      const std::string &outputDir);

  void optimize();
  double sprSearch(unsigned int radius);
  double transferSearch();
  void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);
private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  MiniBMEEvaluator _evaluator;
  std::string _outputDir;
  SpeciesSearchState _searchState;
  bool _rootedScore;
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



