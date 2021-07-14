#pragma once

#include <corax/corax.h>
#include <IO/LibpllParsers.hpp>
#include <string>
#include <memory>
#include <vector>
#include <set>
#include <unordered_set>
#include <util/CArrayRange.hpp>
#include <util/enums.hpp>
#include <util/types.hpp>


typedef struct pll_rnode_s
{
  char * label;
  double length;
  unsigned int node_index;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  struct pll_rnode_s * left;
  struct pll_rnode_s * right;
  struct pll_rnode_s * parent;

  void * data;
} pll_rnode_t;

typedef struct pll_rtree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  pll_rnode_t ** nodes;

  pll_rnode_t * root;

} pll_rtree_t;

void pll_rtree_destroy(pll_rtree_t * tree,
                                  void (*cb_destroy)(void *));


pll_utree_t * pll_rtree_unroot(pll_rtree_t * tree);

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
  pll_rnode_t *getParent(unsigned int node_index) const;
  pll_rnode_t *getNeighbor(unsigned int node_index) const;

  /**
   * labels
   */
  std::unordered_set<std::string> getLabels(bool leavesOnly) const;
  
  static void getLeafLabelsUnder(pll_rnode_t *node,
      std::unordered_set<std::string> &labels);

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

  void ensureUniqueLabels();

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

  /**
   * Get lowest common ancestor
   * First call is O(n^2), and all next calls O(1)
   */
  pll_rnode_t *getLCA(pll_rnode_t *n1, pll_rnode_t *n2);

  /**
   * Return true if either one of n1 or n2 is parent of another
   * First call is O(n^2), and all next calls O(1)
   */
  bool areParents(pll_rnode_t *n1, pll_rnode_t *n2);

  void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);
  
  
  std::vector<bool> &getParentsCache(pll_rnode_t *n1);
  std::vector<bool> &getAncestorssCache(pll_rnode_t *n1);
  
  void buildLCACache();

  /**
   *  Compute and return a mapping between labels and 
   *  unique IDs, such that the mapping does
   *  not depend on the internal indices (to be sure
   *  that the IDs are the same for different processes
   *  working together)
   */
  StringToUint getDeterministicLabelToId() const;
  std::vector<std::string> getDeterministicIdToLabel() const;

  /**
   *  Return a mapping from the nodes IDs of this tree
   *  and the input tree. Both trees must have the same
   *  leaf labels and topology.
   */
  std::vector<unsigned int> getNodeIndexMapping(PLLRootedTree &otherTree);
private:
  std::unique_ptr<pll_rtree_t, void(*)(pll_rtree_t*)> _tree;
  
  struct LCACache {
    // vectors are indexed with rnodes indices
    // lcas[n1][n2] == lca(n1, n2) in the tree
    // parents[n1][n2] is true if n2 is an ancestor
    //   of n1 or n1 an ancestor of n2
    // ancestors[n1][n2] is true if n2 is ancestor of n1
    std::vector<std::vector<pll_rnode_t *> > lcas;
    std::vector<std::vector<bool> > parents;
    std::vector<std::vector<bool> > ancestors;
  };
  std::unique_ptr<LCACache> _lcaCache;
  
  
  static pll_rtree_t *buildRandomTree(const std::unordered_set<std::string> &leafLabels);
};





