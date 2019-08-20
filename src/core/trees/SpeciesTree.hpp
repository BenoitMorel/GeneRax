#pragma once

#include <likelihoods/ReconciliationEvaluation.hpp>
#include <IO/LibpllParsers.hpp>
#include <families/Families.hpp>
#include <string>
#include <maths/DTLRates.hpp>
#include <util/enums.hpp>
#include <memory>

class PerCoreGeneTrees;




class SpeciesTree {
public:
  SpeciesTree(const std::string &newick, bool isFile = true);
  SpeciesTree(const std::unordered_set<std::string> &leafLabels);
  SpeciesTree(const Families &families);
  ~SpeciesTree();
  std::shared_ptr<SpeciesTree> buildRandomTree() const;

  void setRates(const DTLRates &rates);
  void setRatesVector(const DTLRatesVector &rates);
  const DTLRatesVector &getRates() const;
  double computeReconciliationLikelihood(PerCoreGeneTrees &geneTrees, RecModel model);

  std::string toString() const;
  void setRoot(pll_rnode_t *root) {_speciesTree->root = root; root->parent = 0;}
  const pll_rnode_t *getRoot() const {return _speciesTree->root;}
  pll_rnode_t *getRoot() {return _speciesTree->root;}
  unsigned int getTaxaNumber() const;
  pll_rtree_t *getTree() {return _speciesTree;}
  pll_rnode_t *getRandomNode();
  pll_rnode_t *getNode(unsigned int nodeIndex) {return _speciesTree->nodes[nodeIndex];}
  const pll_rnode_t *getNode(unsigned int nodeIndex) const {return _speciesTree->nodes[nodeIndex];}
  unsigned int getMaxNodeIndex() const { return _speciesTree->tip_count + _speciesTree->inner_count;}
  friend std::ostream& operator<<(std::ostream& os, const SpeciesTree &speciesTree) {
    os << speciesTree.toString() << "(" << speciesTree.getTaxaNumber() << " taxa)" << std::endl;
    return os;
  }

  void saveToFile(const std::string &newick, bool masterRankOnly);
private:
  pll_rtree_t *_speciesTree;
  DTLRatesVector _rates;
  void buildFromLabels(const std::unordered_set<std::string> &leafLabels);
  void getLabels(std::unordered_set<std::string> &leafLabels) const;
};


class SpeciesTreeOperator {
public:
  static bool canChangeRoot(const SpeciesTree &speciesTree, int direction);

  /**
   * Change the root to the neighboring branch described by direction where direction is in [0:4[
   */
  static void changeRoot(SpeciesTree &speciesTree, int direction);
  static void revertChangeRoot(SpeciesTree &speciesTree, int direction);
  static unsigned int applySPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft);
  static void reverseSPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int applySPRMoveReturnValue);
  static void getPossiblePrunes(SpeciesTree &speciesTree, std::vector<unsigned int> &prunes);
  static void getPossibleRegrafts(SpeciesTree &speciesTree, unsigned int prune, unsigned int radius, std::vector<unsigned int> &regrafts);
};




