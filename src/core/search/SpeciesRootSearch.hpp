#pragma once

#include <vector>

#include <likelihoods/ReconciliationEvaluation.hpp>
#include <util/types.hpp>

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

class SpeciesTreeLikelihoodEvaluator {
public:
  virtual ~SpeciesTreeLikelihoodEvaluator() {};
  /**
   *  Compute the likelihood of the current species tree
   */
  virtual double computeLikelihood() = 0;

  /**
   *  In case we compute the likelihood on rooted trees,
   *  force gene tree root optimization for the next
   *  computeLikelihood call
   */
  virtual void forceGeneRootOptimization() = 0;

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

  /**
   * Fill perFamLL with the per-family (over all parallel
   * ranks) log-likelihoods
   */
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL) = 0;
};




class SpeciesRootSearch {
public:
  /**
   *  Search for the ML root for the current
   *  species tree topology
   *  
   *  rootLikelihoods and treePerFamLLVec 
   *  will only be filled if not NULL
   */
  static double rootSearch(
      SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluator &evaluation,
      unsigned int maxDepth,
      RootLikelihoods *rootLikelihoods = nullptr,
      TreePerFamLLVec *treePerFamLLVec = nullptr); 
};

