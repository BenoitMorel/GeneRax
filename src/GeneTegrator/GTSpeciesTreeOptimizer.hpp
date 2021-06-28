#pragma once

#include <search/SpeciesRootSearch.hpp>
#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"
#include <trees/PLLRootedTree.hpp>
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
  void setEvaluations(PerCoreMultiEvaluation &evaluations) {
    _evaluations = &evaluations;
  }
  virtual ~GTSpeciesTreeLikelihoodEvaluator() {}
  virtual double computeLikelihood();
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const {return false;}
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL);
  virtual void countEvents();
  virtual void getTransferInformation(PLLRootedTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents);
  virtual bool pruneSpeciesTree() const {assert(false); return false;}
private:
  PerCoreMultiEvaluation *_evaluations;
};


class GTSpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  GTSpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &families, 
      const RecModelInfo &info,
      const std::string &outputDir);

  double computeRecLikelihood();
  double sprSearch(unsigned int radius);
  double rootSearch(unsigned int maxDepth);

  void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);

private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  PerCoreMultiEvaluation _evaluations;
  GTSpeciesTreeLikelihoodEvaluator _evaluator;
  std::string _outputDir;
  double _bestRecLL;

  bool testPruning(unsigned int prune,
    unsigned int regraft);
  void newBestTreeCallback(double newLL);
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



