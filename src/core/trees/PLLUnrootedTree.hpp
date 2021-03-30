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

#include <unordered_set>
#include <string>
#include <memory>
#include <vector>
#include <unordered_set>
#include <util/types.hpp>
#include <util/CArrayRange.hpp>
#include <maths/Random.hpp>
#include <functional>

class PLLRootedTree;


using UnodePrinter = 
  std::function<void(pll_unode_t *, std::stringstream &)>;

void defaultUnodePrinter(pll_unode_t *node, 
    std::stringstream &ss);

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


  PLLUnrootedTree(PLLRootedTree &rootedTree);

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
  
  std::unordered_set<std::string> getLabels() const;

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
  std::string getNewickString(UnodePrinter f = defaultUnodePrinter,
      pll_unode_t *root = nullptr, 
      bool rooted = false);

  /**
   * Replace null branch lengths with minBL
   */
  void setMissingBranchLengths(double minBL = 0.1); 

  /**
   *  Direct access to the libpll structure
   */
  pll_utree_t *getRawPtr() {return _tree.get();}
  
  CArrayRange<pll_unode_t*> getLeaves() const;

  /*
   *  C++11 range for accessing internal nodes. 
   *  Only includes one of the three pll elements per node
   *  in the tree
   */
  CArrayRange<pll_unode_t*> getInnerNodes();

  /*
   *  C++11 range for accessing nodes. 
   *  For internal nodes, Only includes one of the three pll elements per node
   *  in the tree
   */
  CArrayRange<pll_unode_t*> getNodes();

  /**
   *  Create a vector of all nodes (including all three pll
   *  internal elements per internal node), such that a node
   *  always comes after its (virtual) children.
   */
  std::vector<pll_unode_t*> getPostOrderNodes(bool innerOnly = false);

  /**
   *  Return a set of all branches (for each node, one and only 
   *  one of node and node->back will be inserted in the set)
   */
  std::unordered_set<pll_unode_t *> getBranches();

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


  /**
   *  replace pu and pv by one of their three possible next pointers
   *  such that their back pointers point to each other
   *  In addition, fill branchesPath with the sequence of BRANCHES
   *  along the path (if node1->back == node2, we only add one of the two)
   */
  static void orientTowardEachOther(pll_unode_t **pu, 
      pll_unode_t **pv,
      std::vector<pll_unode_t *> &branchesPath);


  /*
   *  Return the virtual root in the unrooted tree corresponding to
   *  the root in the input rooted tree. Both trees should have the 
   *  same (unrooted) topology and leaf set.
   */
   pll_unode_t *getVirtualRoot(PLLRootedTree &referenceTree);

  // TODO: DOCUMENT
  std::vector<double> getMADRelativeDeviations();

  /**
   *  Return the set of leaves under the input directed node
   */
  static std::unordered_set<unsigned int> getClade(pll_unode_t *node);

private:
  std::unique_ptr<pll_utree_t, void(*)(pll_utree_t*)> _tree;
};




