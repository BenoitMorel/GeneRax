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


typedef struct corax_rnode_s
{
  char * label;
  double length;
  unsigned int node_index;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  struct corax_rnode_s * left;
  struct corax_rnode_s * right;
  struct corax_rnode_s * parent;

  void * data;
} corax_rnode_t;

typedef struct corax_rtree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;

  corax_rnode_t ** nodes;

  corax_rnode_t * root;

} corax_rtree_t;

void corax_rtree_destroy(corax_rtree_t * tree,
                                  void (*cb_destroy)(void *));


corax_utree_t * corax_rtree_unroot(corax_rtree_t * tree);

/**
 *  C++ wrapper around the libpll corax_rtree_t structure
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
  corax_rnode_t *getRoot() const;
  corax_rnode_t *getAnyInnerNode() const;
  corax_rnode_t *getNode(unsigned int node_index) const;
  corax_rnode_t *getParent(unsigned int node_index) const;
  corax_rnode_t *getNeighbor(unsigned int node_index) const;

  /**
   * labels
   */
  std::unordered_set<std::string> getLabels(bool leavesOnly) const;
  
  static void getLeafLabelsUnder(corax_rnode_t *node,
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
  corax_rtree_t *getRawPtr() {return _tree.get();}
  const corax_rtree_t *getRawPtr() const {return _tree.get();}
  
  CArrayRange<corax_rnode_t*> getLeaves() const;
  CArrayRange<corax_rnode_t*> getInnerNodes() const;
  CArrayRange<corax_rnode_t*> getNodes() const;
  std::vector<corax_rnode_t*> getPostOrderNodes() const;

  static void setSon(corax_rnode_t *parent, corax_rnode_t *newSon, bool left);
  
  friend std::ostream& operator<<(std::ostream& os, const PLLRootedTree &tree)
  {
    char *newick = corax_rtree_export_newick(tree.getRawPtr()->root, 0);
    std::string str(newick);
    os << str;
    free(newick);
    return os;
  }

  /**
   * Get lowest common ancestor
   * First call is O(n^2), and all next calls O(1)
   */
  corax_rnode_t *getLCA(corax_rnode_t *n1, corax_rnode_t *n2);
  corax_rnode_t *getLCA(unsigned int nodeIndex1, unsigned int nodeIndex2);

  /**
   * Return true if either one of n1 or n2 is parent of another
   * First call is O(n^2), and all next calls O(1)
   */
  bool areParents(corax_rnode_t *n1, corax_rnode_t *n2);

  void onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  
  
  std::vector<bool> &getParentsCache(corax_rnode_t *n1);
  std::vector<bool> &getAncestorssCache(corax_rnode_t *n1);
  
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

  bool areNodeIndicesParallelConsistent() const;
private:
  std::unique_ptr<corax_rtree_t, void(*)(corax_rtree_t*)> _tree;
  
  struct LCACache {
    // vectors are indexed with rnodes indices
    // lcas[n1][n2] == lca(n1, n2) in the tree
    // parents[n1][n2] is true if n2 is an ancestor
    //   of n1 or n1 an ancestor of n2
    // ancestors[n1][n2] is true if n2 is ancestor of n1
    std::vector<std::vector<corax_rnode_t *> > lcas;
    std::vector<std::vector<bool> > parents;
    std::vector<std::vector<bool> > ancestors;
  };
  std::unique_ptr<LCACache> _lcaCache;
  
  
  static corax_rtree_t *buildRandomTree(const std::unordered_set<std::string> &leafLabels);
};





