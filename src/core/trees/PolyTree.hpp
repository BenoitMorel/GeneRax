#pragma once

#include <vector>
#include <string>

#include <trees/PLLUnrootedTree.hpp>

class PolyTree
{
  public:
    PolyTree(PLLUnrootedTree &tree);

    std::string getNewickString() const;

    const std::vector<pll_unode_t *> & getChildren(pll_unode_t *node) const {
      return _cells[node->node_index].children;;
    }

    struct Cell {
      unsigned int gid;
      std::vector<pll_unode_t *> children;
    };
  private:
    PLLUnrootedTree &_tree;
    std::vector<Cell> _cells; 
};


