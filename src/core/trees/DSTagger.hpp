#pragma once

#include <vector>
#include <set>
#include <trees/PLLUnrootedTree.hpp>
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
  bool goesUp(unsigned int nodeIndex) const {
    return _clvs[nodeIndex].goesUp;
  }
  void print();
private:
  PLLUnrootedTree &_tree;
  
  using Clade = std::set<unsigned int>;
  struct CLV {
    unsigned int score;
    bool isDup;
    Clade clade;
    bool goesUp;
    CLV(): score(0), isDup(false), goesUp(false) {}
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



