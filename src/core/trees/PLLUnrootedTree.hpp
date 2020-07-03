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
#include <util/types.hpp>
#include <util/CArrayRange.hpp>
#include <maths/Random.hpp>

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
  PLLUnrootedTree(const std::vector<const char*> &labels, unsigned int seed = static_cast<unsigned int>(Random::getInt()));

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
  // number of nodes in the node buffer (#leaves + #innner)
  unsigned int getNodesNumber() const;
  // #leaves + 3 * #inner
  unsigned int getDirectedNodesNumber() const;
  // #leaves
  unsigned int getLeavesNumber() const;
  // #inner
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

  /*
   *  C++11 range for accessing nodes. 
   *  Only includes one of the three pll internal element per node
   *  in the tree
   */
  CArrayRange<pll_unode_t*> getNodes();

  /**
   *  Create a vector of all nodes (including all three pll
   *  internal elements per internal node), such that a node
   *  always comes after its (virtual) children.
   */
  std::vector<pll_unode_t*> getPostOrderNodes();

  /**
   *  Compute a matrix of pairwise distances.
   *  The distance between two nodes is the sum of the
   *  branch lengths along the path between these two nodes.
   *  
   *  distances[i][j] is the distance between node1 and node2,
   *  such that node1->node_index == i && node2->node_index == j
   *
   *  If leavesOnly is true, the output matrix is symetric, 
   *  and the diagonal elements are set to 0.0
   *
   *  If leavesOnly is false, distances is a rectangular matrix
   *  with dimensions distances[max_gene_id][leaves_number]
   *
   *  
   */
  void computePairwiseDistances(MatrixDouble &distances,
      bool leavesOnly = true);

  // TODO: DOCUMENT
  std::vector<double> getMADRelativeDeviations();

  /**
   *  Return the set of leaves under the input directed node
   */
  static std::unordered_set<unsigned int> getClade(pll_unode_t *node);

private:
  std::unique_ptr<pll_utree_t, void(*)(pll_utree_t*)> _tree;
};




