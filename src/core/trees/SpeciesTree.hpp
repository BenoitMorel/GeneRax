#pragma once

#include <likelihoods/ReconciliationEvaluation.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/Families.hpp>
#include <string>
#include <util/enums.hpp>
#include <memory>
#include <trees/PLLRootedTree.hpp>
#include <unordered_set>
#include <unordered_map>

class PerCoreGeneTrees;

class SpeciesTree {
public:
  SpeciesTree(const std::string &newick, bool isFile = true);
  SpeciesTree(const std::unordered_set<std::string> &leafLabels);
  SpeciesTree(const Families &families);
  // forbid copy
  SpeciesTree(const SpeciesTree &) = delete;
  SpeciesTree & operator = (const SpeciesTree &) = delete;
  SpeciesTree(SpeciesTree &&) = delete;
  SpeciesTree & operator = (SpeciesTree &&) = delete;
  
  std::unique_ptr<SpeciesTree> buildRandomTree() const;

  // set rates and gene trees before calling this

  std::string toString() const;
  pll_rnode_t *getRandomNode();
  pll_rnode_t *getNode(unsigned int nodeIndex) {return _speciesTree.getNode(nodeIndex);}
  pll_rnode_t *getRoot() { return getTree().getRawPtr()->root; } 
  friend std::ostream& operator<<(std::ostream& os, SpeciesTree &speciesTree) {
    os << speciesTree.toString() << "(" << speciesTree.getTree().getLeavesNumber() << " taxa)" << std::endl;
    return os;
  }

  const PLLRootedTree &getTree() const {return _speciesTree;}
  PLLRootedTree &getTree() {return _speciesTree;}

  void saveToFile(const std::string &newick, bool masterRankOnly);
  size_t getHash() const;
  size_t getNodeIndexHash() const;
  void getLabelsToId(std::unordered_map<std::string, unsigned int> &map) const;

  class Listener {
  public:
    virtual ~Listener() {}
    virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate) = 0;
  };
  void addListener(Listener *listener);
  void removeListener(Listener *listener);
  void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate); // should be called when changing the species tree


private:
  PLLRootedTree _speciesTree;
  std::vector<Listener *> _listeners;
  void buildFromLabels(const std::unordered_set<std::string> &leafLabels);
  static std::unordered_set<std::string> getLabelsFromFamilies(const Families &families);
};


class SpeciesTreeOperator {
public:
  static bool canChangeRoot(const SpeciesTree &speciesTree, unsigned int direction);

  /**
   * Change the root to the neighboring branch described by direction where direction is in [0:4[
   */
  static void changeRoot(SpeciesTree &speciesTree, unsigned int direction);
  static void revertChangeRoot(SpeciesTree &speciesTree, unsigned int direction);
  static bool canApplySPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft);
  static unsigned int applySPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft);
  static void reverseSPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int applySPRMoveReturnValue);
  static void getPossiblePrunes(SpeciesTree &speciesTree, 
      std::vector<unsigned int> &prunes, 
      std::vector<double> support, 
      double maxSupport);
  static void getPossibleRegrafts(SpeciesTree &speciesTree, 
      unsigned int prune, 
      unsigned int radius, 
      std::vector<unsigned int> &regrafts);

};




