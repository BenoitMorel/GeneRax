#pragma once

#include <IO/LibpllParsers.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <maths/DTLRates.hpp>
#include <util/enums.hpp>

class PerCoreGeneTrees;




class SpeciesTree {
public:
  SpeciesTree(const std::string &newick, bool isFile = true); 
  ~SpeciesTree();

  void setRates(const DTLRates &rates);
  const DTLRates &getRates() const;
  double computeReconciliationLikelihood(PerCoreGeneTrees &geneTrees, RecModel model);


  std::string toString() const;

  friend std::ostream& operator<<(std::ostream& os, const SpeciesTree &speciesTree) {
    os << speciesTree.toString() << "(" << speciesTree.getTaxaNumber() << " taxa)" << std::endl;
    return os;
  }

  void setRoot(pll_rnode_t *root) {_speciesTree->root = root; root->parent = 0;}
  const pll_rnode_t *getRoot() const {return _speciesTree->root;}
  pll_rnode_t *getRoot() {return _speciesTree->root;}
  unsigned int getTaxaNumber() const;

private:
  pll_rtree_t *_speciesTree;
  DTLRates _rates;
};


class SpeciesTreeOperator {
public:
  static bool canChangeRoot(const SpeciesTree &speciesTree, bool left1);

  /**
   * Change the root to the branch between root->side1 and root->side1->side2
   * where side can be left or right
   * Can be reverted with !left1 and !left2
   */
  static void changeRoot(SpeciesTree &speciesTree, bool left1, bool left2);
  
};

class SpeciesTreeOptimizer {
};

