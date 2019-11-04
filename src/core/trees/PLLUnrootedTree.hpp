#pragma once

extern "C" {
#include <pll.h>
#include <pllmod_algorithm.h>
#include <pll_binary.h>
#include <pll_msa.h>
#include <pll_optimize.h>
#include <pll_tree.h>
#include <pllmod_util.h>  
}

#include <string>
#include <memory>
#include <vector>
#include <unordered_set>
#include <util/CArrayRange.hpp>


/**
 *  C++ wrapper around the libpll pll_utree_t structure
 *  to represent an unrooted tree
 */
class PLLUnrootedTree {
public:
  /**
   * Construct from a string that is either a path 
   * to a newick file or a newick string
   */
  PLLUnrootedTree(const std::string &str, bool isFile = true);

  /**
   *  Construct a random tree from a set of taxa labels
   */
  PLLUnrootedTree(const std::vector<const char*> &labels, unsigned int seed = rand());

  /**
   * Forbid copy
   */
  PLLUnrootedTree(const PLLUnrootedTree &) = delete;
  PLLUnrootedTree & operator = (const PLLUnrootedTree &) = delete;
  PLLUnrootedTree(PLLUnrootedTree &&) = delete;
  PLLUnrootedTree & operator = (PLLUnrootedTree &&) = delete;


  /*
   * Tree dimension
   */
  unsigned int getNodesNumber() const;
  unsigned int getLeavesNumber() const;
  unsigned int getInnerNodesNumber() const; 

  /*
   * Node access
   */
  pll_unode_t *getAnyInnerNode();
  pll_unode_t *getNode(unsigned int node_index);

 
  
  /**
   * labels
   */
  std::unordered_set<std::string> getLeavesLabels();

  /*
   * Save the tree in newick format in filename
   */
  void save(const std::string &fileName); 

  /**
   * Replace null branch lengths with minBL
   */
  void setMissingBranchLengths(double minBL = 0.1); 

  /**
   *  Direct access to the libpll structure
   */
  pll_utree_t *getRawPtr() {return _tree.get();}
  
  CArrayRange<pll_unode_t*> getLeaves();
private:
  std::unique_ptr<pll_utree_t, void(*)(pll_utree_t*)> _tree;
};




