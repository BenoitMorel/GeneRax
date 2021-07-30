#pragma once

#include <ccp/SpeciesSplits.hpp>
#include <search/SpeciesRootSearch.hpp>
#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>
#include <vector>
#include <maths/ModelParameters.hpp>
#include <ccp/SpeciesSplitScore.hpp>


class SpeciesSplitEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  SpeciesSplitEvaluator(PLLRootedTree &rootedTree,
      const Families &families,
      bool rootedScore);
  virtual ~SpeciesSplitEvaluator() {}
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
  SpeciesSplits _splits;
  std::unique_ptr<SpeciesSplitScore> _score;
};


class SpeciesSplitSpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  SpeciesSplitSpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &families, 
      const std::string &outputDir,
      bool rootedScore);

  void optimize();
  double sprSearch(unsigned int radius);
  double transferSearch();
  void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);
  double rootSearch(unsigned int maxDepth);
private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  SpeciesSplitEvaluator _evaluator;
  std::string _outputDir;
  SpeciesSearchState _searchState;
  bool _rootedScore;
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



