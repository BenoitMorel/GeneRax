#pragma once

#include <vector>

#include <likelihoods/ReconciliationEvaluation.hpp>
#include <util/types.hpp>
#include <util/Scenario.hpp>
#include <maths/AverageStream.hpp>

class SpeciesTree; 

using TreePerFamLL = std::pair<std::string, PerFamLL>;
using TreePerFamLLVec = std::vector<TreePerFamLL>;
struct RootLikelihoods {
  void reset() {
    idToLL.clear();
  }
  void saveValue(pll_rnode_t *t, double ll);
  void fillTree(PLLRootedTree &tree);
  std::unordered_map<std::string, double> idToLL;
};

/**
 *  Interface for classes that are used by the search
 *  algorithms to evaluate the reconciliation likelihood and 
 *  to describe the reconciliation model.
 */
class SpeciesTreeLikelihoodEvaluatorInterface {
public:
  virtual ~SpeciesTreeLikelihoodEvaluatorInterface() {};
  /**
   *  Exact likelihood computation
   */
  virtual double computeLikelihood() = 0;

  /**
   *  Fast but approximated version of the likelihood
   *  computation
   */
  virtual double computeLikelihoodFast() = 0;

  /**
   *  Returns true if computeLikelihood and computeLikelihoodFast
   *  have different implementations
   */
  virtual bool providesFastLikelihoodImpl() const = 0;


  /**
   *  Save in a stack what needs to be saved in case 
   *  of the rollback of a species tree operation
   */
  virtual void pushRollback() = 0;

  /**
   * Pop the upper state and apply it after a rollback
   * of a species tree operation
   */
  virtual void popAndApplyRollback() = 0;
  
  virtual void getTransferInformation(PLLRootedTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents) = 0;
  
  /**
   * Fill perFamLL with the per-family (over all parallel
   * ranks) log-likelihoods
   */
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL) = 0;
  
  /**
   *  Are we in prune species tree mode?
   */
  virtual bool pruneSpeciesTree() const = 0;
};

/**
 *  Structure that describes the current state of the
 *  species tree search, to decide which strategy to adopt:
 *  for instance, when the tree is far from being plausible,
 *  we apply a less thorough but faster search strategy.
 */
struct SpeciesSearchState {
public:
  SpeciesSearchState(): farFromPlausible(true) {}
  /**
   *  Set to true when the tree is far from being plausible
   *  (in the early stage of the search after starting from
   *  a random tree for instance).
   *  
   *  When set to true, the search strategy can choose to
   *  apply less thorough optimizations and to favor speed
   *  over small tree improvments.
   *
   *  Both the search algorithm functions and their caller 
   *  can update this variable.
   */
  bool farFromPlausible;

  /**
   *  Stores the average difference between the average and
   *  the exact likelihood values. This average is updated
   *  in a streaming fashion everytime that both the average
   *  and exact likelihood are computed during the search
   *  (see for instance SpeciesSearchCommon::testSPR)
   *  It is then used to decide if the approximated likelihood
   *  score is good enough to try estimating the exact score.
   *  
   *  This is only relevant when 
   *  SpeciesTreeLikelihoodEvaluatorInterface::providesFastLikelihoodImpl
   *  is set to true
   */
  AverageStream averageApproxError;
};

class SpeciesSearchCommon {
public:
  /**
   *  Test a SPR move. 
   *  If it improves the likelihood, keep it and return true.
   *  Else, rollback it and return false.
   */
  static bool testSPR(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    unsigned int prune,
    unsigned int regraft,
    SpeciesSearchState &searchState,
    double bestLL,
    double &testedTreeLL
    );

  /**
   *  Try SPR moves with a small radius around the species
   *  node with id spid. If a move improves the likelihood, 
   *  apply it, and recursively search further. Otherwise, 
   *  rollback the tested moves.
   *  Returns true if one better tree has been found
   */
  static bool veryLocalSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int spid,
    double previousBestLL,
    double &newBestLL);
  
};

