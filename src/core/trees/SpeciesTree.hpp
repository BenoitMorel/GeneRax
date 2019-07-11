#pragma once

#include <IO/LibpllParsers.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <maths/DTLRates.hpp>
#include <util/enums.hpp>

class PerCoreGeneTrees;




class SpeciesTree {
public:
  SpeciesTree(const std::string &newick, bool isFile = true);
  SpeciesTree(const std::unordered_set<std::string> &leafLabels);
  SpeciesTree(const std::vector<FamiliesFileParser::FamilyInfo> &families);
  ~SpeciesTree();

  void setRates(const DTLRates &rates);
  const DTLRates &getRates() const;
  double computeReconciliationLikelihood(PerCoreGeneTrees &geneTrees, RecModel model);

  std::string toString() const;
  void setRoot(pll_rnode_t *root) {_speciesTree->root = root; root->parent = 0;}
  const pll_rnode_t *getRoot() const {return _speciesTree->root;}
  pll_rnode_t *getRoot() {return _speciesTree->root;}
  unsigned int getTaxaNumber() const;
  
  friend std::ostream& operator<<(std::ostream& os, const SpeciesTree &speciesTree) {
    os << speciesTree.toString() << "(" << speciesTree.getTaxaNumber() << " taxa)" << std::endl;
    return os;
  }

  void saveToFile(const std::string &newick);
private:
  pll_rtree_t *_speciesTree;
  DTLRates _rates;
  void buildFromLabels(const std::unordered_set<std::string> &leafLabels);
};


class SpeciesTreeOperator {
public:
  static bool canChangeRoot(const SpeciesTree &speciesTree, int direction);

  /**
   * Change the root to the neighboring branch described by direction where direction is in [0:4[
   */
  static void changeRoot(SpeciesTree &speciesTree, int direction);
  static void revertChangeRoot(SpeciesTree &speciesTree, int direction);
};




class SpeciesTreeOptimizer {
public:
  static void rootSlidingSearch(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model);
  static void rootExhaustiveSearch(SpeciesTree &speciesTree, PerCoreGeneTrees &geneTrees, RecModel model);
};
