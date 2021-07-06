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

class RecModelInfo;
using MultiEvaluation = ReconciliationModelInterface;
using MultiEvaluationPtr = 
  std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluation = std::vector<MultiEvaluationPtr>;

class GTSpeciesTreeLikelihoodEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  GTSpeciesTreeLikelihoodEvaluator()
  {}
  void setEvaluations(const RecModelInfo &info, 
      const Families &families,
      PerCoreMultiEvaluation &evaluations,
      PerCoreGeneTrees &geneTrees) {
    _info = info;
    _families = &families;
    _evaluations = &evaluations;
    _geneTrees = &geneTrees;
  }
  virtual ~GTSpeciesTreeLikelihoodEvaluator() {}
  virtual double computeLikelihood();
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const {return false;}
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL);
  virtual void getTransferInformation(PLLRootedTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents);
  virtual bool pruneSpeciesTree() const {return _info.pruneSpeciesTree;}
private:
  RecModelInfo _info;
  const Families *_families;
  PerCoreMultiEvaluation *_evaluations;
  PerCoreGeneTrees *_geneTrees;
};


class GTSpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  GTSpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &families, 
      const RecModelInfo &info,
      const std::string &outputDir);

  void optimize();
  double sprSearch(unsigned int radius);
  double rootSearch(unsigned int maxDepth);
  double transferSearch();
  void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);

private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  PerCoreGeneTrees _geneTrees;
  PerCoreMultiEvaluation _evaluations;
  GTSpeciesTreeLikelihoodEvaluator _evaluator;
  RecModelInfo _info;
  std::string _outputDir;
  SpeciesSearchState _searchState;
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



