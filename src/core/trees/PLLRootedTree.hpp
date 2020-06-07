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
#include <util/enums.hpp>


/**
 *  C++ wrapper around the libpll pll_rtree_t structure
 *  to represent a rooted tree
 */
class PLLRootedTree {
public:
  /**
   * Construct from a string that is either a path 
   * to a newick file or a newick string
   */
  PLLRootedTree(const std::string &str, bool isFile = true);

  /**
   *  Construct a random tree from a set of taxa labels
   */
  PLLRootedTree(const std::unordered_set<std::string> &labels);

  /**
   * Forbid copy
   */
  PLLRootedTree(const PLLRootedTree &) = delete;
  PLLRootedTree & operator = (const PLLRootedTree &) = delete;
  PLLRootedTree(PLLRootedTree &&) = delete;
  PLLRootedTree & operator = (PLLRootedTree &&) = delete;


  /*
   * Tree dimension
   */
  unsigned int getNodesNumber() const;
  unsigned int getLeavesNumber() const;
  unsigned int getInnerNodesNumber() const; 

  /*
   * Node access
   */
  pll_rnode_t *getRoot() const;
  pll_rnode_t *getAnyInnerNode() const;
  pll_rnode_t *getNode(unsigned int node_index) const;

  /**
   * labels
   */
  std::unordered_set<std::string> getLabels(bool leavesOnly) const;

  /**
   *  Get a mapping from label to integer 
   */
  StringToUintMap getLabelToIntMap();

  /*
   * Save the tree in newick format in filename
   */
  void save(const std::string &fileName) const; 

  std::string getNewickString() const;

  /**
   * Replace null branch lengths with minBL
   */
  void setMissingBranchLengths(double minBL = 1.0); 

  void setMissingLabels();

  /**
   *  Direct access to the libpll structure
   */
  pll_rtree_t *getRawPtr() {return _tree.get();}
  const pll_rtree_t *getRawPtr() const {return _tree.get();}
  
  CArrayRange<pll_rnode_t*> getLeaves() const;
  CArrayRange<pll_rnode_t*> getInnerNodes() const;
  CArrayRange<pll_rnode_t*> getNodes() const;
  std::vector<pll_rnode_t*> getPostOrderNodes() const;

  static void setSon(pll_rnode_t *parent, pll_rnode_t *newSon, bool left);
  
  friend std::ostream& operator<<(std::ostream& os, const PLLRootedTree &tree)
  {
    char *newick = pll_rtree_export_newick(tree.getRawPtr()->root, 0);
    std::string str(newick);
    os << str;
    free(newick);
    return os;
  }



private:
  std::unique_ptr<pll_rtree_t, void(*)(pll_rtree_t*)> _tree;

  static pll_rtree_t *buildRandomTree(const std::unordered_set<std::string> &leafLabels);
};





