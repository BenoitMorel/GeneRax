#pragma once

#include <search/SpeciesRootSearch.hpp>
#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"
#include <trees/PLLRootedTree.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <memory>
#include <vector>
#include <maths/ModelParameters.hpp>

class RecModelInfo;
using MultiEvaluation = MultiModel;
using MultiEvaluationPtr = 
  std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluation = std::vector<MultiEvaluationPtr>;

class GTSpeciesTreeLikelihoodEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  GTSpeciesTreeLikelihoodEvaluator(SpeciesTree &speciesTree,
      ModelParameters &modelRates, 
      const Families &families,
      PerCoreGeneTrees &geneTrees);
  virtual ~GTSpeciesTreeLikelihoodEvaluator() {}
  virtual double computeLikelihood();
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const {return false;}
  virtual bool isDated() const {return _modelRates.info.isDated();}
  virtual double optimizeModelRates(bool thorough = false);
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL);
  virtual void getTransferInformation(SpeciesTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents);
  virtual bool pruneSpeciesTree() const {return _modelRates.info.pruneSpeciesTree;}
  virtual void setAlpha(double alpha);
  
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  void printHightPrecisionCount();
protected:
  virtual double optimizeGammaRates();
private:
  SpeciesTree &_speciesTree;
  ModelParameters &_modelRates;
  const Families &_families;
  PerCoreMultiEvaluation _evaluations;
  PerCoreMultiEvaluation _approxEvaluations;
  std::vector<int> _highPrecisions;
  PerCoreGeneTrees &_geneTrees;
  
};


class GTSpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  GTSpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &families, 
      const RecModelInfo &info,
      const std::string &outputDir);

  void optimize();
  double sprSearch(unsigned int radius);
  double rootSearch(unsigned int maxDepth, bool thorough = false);
  double transferSearch();
  void onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  void reconcile(unsigned int samples);
  double optimizeModelRates(bool thorough = false);
  void optimizeDates(bool thorough = true);
  SpeciesTree &getSpeciesTree() {return *_speciesTree;}
  void randomizeRoot();
  void saveSpeciesTree();
  void saveSpeciesTreeRootLL();
  SpeciesTreeLikelihoodEvaluatorInterface &getEvaluator() {return *_evaluator;}
private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  PerCoreGeneTrees _geneTrees;
  RecModelInfo _info;
  PerCoreMultiEvaluation _evaluations;
  std::unique_ptr<GTSpeciesTreeLikelihoodEvaluator> _evaluator;
  ModelParameters _modelRates;
  std::string _outputDir;
  SpeciesSearchState _searchState;
  RootLikelihoods _rootLikelihoods;
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



