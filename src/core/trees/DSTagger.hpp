#pragma once

#include <vector>
#include <set>
#include <trees/PLLUnrootedTree.hpp>
#include <util/types.hpp>
class Clade;

/**
 *  Apply astral-pro tagging
 *  Assumes that the species ID is stored in 
 *  the clv_index field of the pll_unode_t structure
 */
class DSTagger {
public:
  DSTagger(PLLUnrootedTree &tree);
  pll_unode_t * &getRoot() {
    return _bestRoots[0];
  }
  bool isDuplication(unsigned int nodeIndex) const {
    return _clvs[nodeIndex].isDup;
  }
  void print();

  /*
   *  Go up in the rooted tree from node (not included)
   *  At each traversed parent, go down on the side
   *  not leading to node, and fill set with all
   *  subnodes clv_index (we assume that 
   *  clv_index contains the SPID)
   *  We skip parents tagged as duplications
   */
  void fillUpTraversal(pll_unode_t *node,
      TaxaSet &set);

  /**
   *  Fill set with all nodes under node,
   *  where "under" is defined from the root
   *  stored in the DSTagger instance.
   */
  void fillWithChildren(pll_unode_t *node,
      TaxaSet &set);

  /**
   *  Fills descendants with all nodes under node,
   *  node should already be oriented (toward the root) 
   */
  void fillWithInternalDescendants(pll_unode_t *node,
      std::vector<pll_unode_t*> &descendants);

  pll_unode_t *getThirdNode(pll_unode_t *n1,
      pll_unode_t *n2) {
    return (n1->next == n2) ? n2->next : n1->next;
  }


  /**
   * Get all ancestors of node that are not 
   * duplications and that are oriented toward node
   */
  std::vector<pll_unode_t *> getSpeciationAncestorNodes(pll_unode_t *node);

  /**
   *  returns the node (among the three representing the 
   *  same internal node) pointing toward the root
   */
  void orientUp(pll_unode_t *&node) const {
    node = _clvs[node->node_index].up;
  }

  bool goesUp(pll_unode_t *node) const {
    return node == _clvs[node->node_index].up;
  }
private:
  PLLUnrootedTree &_tree;
  bool _isRootDup; 

  std::array<pll_unode_t, 3> _roots;
  using Clade = std::set<unsigned int>;
  struct CLV {
    unsigned int score;
    bool isDup;
    bool isRoot;
    Clade clade;
    pll_unode_t *up;
    CLV(): score(0), 
      isDup(false), 
      isRoot(false), 
      up(nullptr)
      {}
  };
  
  std::vector<CLV> _clvs;
  std::vector<pll_unode_t *> _bestRoots;
  void _tagNode(pll_unode_t *node, CLV &clv);
  void _rootFromNode(pll_unode_t *node);

  struct TaggerUNodePrinter {
    TaggerUNodePrinter(std::vector<CLV> &clvs):clvs(clvs) { }
    std::vector<CLV> &clvs;
    void operator()(pll_unode_t *node, 
        std::stringstream &ss) const;
  };
};



