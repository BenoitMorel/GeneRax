#pragma once

#include <trees/PLLRootedTree.hpp>
#include <vector>

class DatedTree {
public:
  DatedTree(PLLRootedTree *rootedTree, bool fromBL = false);
  


  void rescaleBranchLengths();

  PLLRootedTree &getRootedTree() const {return *_rootedTree;}
  const std::vector<corax_rnode_s *> &getOrderedSpeciations() const {return _orderedSpeciations;}
  unsigned int getRank(unsigned int nodeIndex) const {return _ranks[nodeIndex];}
  const std::vector<unsigned int> &getOrderedSpeciesRanks() const {return _ranks;}


  bool moveDown(unsigned int rank, bool force = false);
  bool moveUp(unsigned int rank, bool force = false);

  void moveNodeToRoot(corax_rnode_t *node);

  bool isConsistent() const;

  struct Backup {
    std::vector<unsigned int> ranks;
  };

  Backup getBackup() {Backup backup; backup.ranks = _ranks; return backup;}
  void restore(const Backup &backup);

  /**
   *  hash value that characterizes the current order of the speciation events
   *
   */
  size_t getOrderingHash(size_t startingHash = 42) const;
private:
  PLLRootedTree *_rootedTree;
 
  // internal nodes, from the root to the most recent speciation 
  std::vector<corax_rnode_s *> _orderedSpeciations;
  // parents have always a lower rank than their children
  std::vector<unsigned int> _ranks;
};


