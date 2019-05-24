#pragma once

#include <unordered_map>
#include <unordered_set>
#include <likelihoods/LibpllEvaluation.hpp>
#include <vector>
#include <IO/Logger.hpp>


class SubtreeRepeatsCache;

struct hashing_func_subtree {
  hashing_func_subtree(SubtreeRepeatsCache &cache): _cache(cache) {}
  unsigned long operator()(const pll_unode_t *subtree) const;
  SubtreeRepeatsCache &_cache;
};

struct key_equal_fn_subtree {
  key_equal_fn_subtree(SubtreeRepeatsCache &cache): _cache(cache) {}
  bool operator()(const pll_unode_t *subtree1, const pll_unode_t *subtree2) const;
  SubtreeRepeatsCache &_cache;
};


static void fillPreOrderRec(pll_unode_t *node, std::vector<pll_unode_t *> &nodes, std::unordered_set<pll_unode_t *> &marked)
{
  if (marked.find(node) != marked.end()) {
    return;
  }
  marked.insert(node);
  if (node->next) {
    fillPreOrderRec(node->next->back, nodes, marked);  
    fillPreOrderRec(node->next->next->back, nodes, marked);  
  }

  nodes.push_back(node);
}

static void fillPreOrder(pll_utree_t *tree, std::vector<pll_unode_t *> &nodes)
{
  std::unordered_set<pll_unode_t *> marked;
  nodes.clear();
  auto nodesNumber = tree->tip_count + tree->inner_count;
  for (unsigned int i = 0; i < nodesNumber; ++i) {
    auto node = tree->nodes[i];
    fillPreOrderRec(node, nodes, marked);
    if (node->next) {
      fillPreOrderRec(node->next, nodes, marked);
      fillPreOrderRec(node->next->next, nodes, marked);
    }
  }
}

class SubtreeRepeatsCache {
public:
  SubtreeRepeatsCache(): 
    _subtreeToRID(10, hashing_func_subtree(*this), key_equal_fn_subtree(*this)), 
    _enabled(false)
  {}

  void setGenesToSpecies(const std::vector<unsigned int> &geneToSpecies) {
    _geneToSpecies = geneToSpecies;
  }

  void resetCache() {
    _subtreeToRID.clear();
    _RIDToSubtree.clear();
    _NIDToRID = std::vector<unsigned int>(_geneToSpecies.size() + 1, static_cast<unsigned int>(-1));
  }

  pll_unode_t *getRepeat(pll_unode_t *subtree) {
    if (!_enabled) {
      return subtree;
    }
    auto repeatIndex = getRepeatIndexNoCheck(subtree);
    return _RIDToSubtree[repeatIndex];
  }
  
  void setTree(pll_utree_t *tree) {
    if (_subtreeToRID.size() > 0 || !_enabled) {
      return;
    }
    resetCache();
    std::vector<pll_unode_t *> preOrderNodes;
    fillPreOrder(tree, preOrderNodes);
    assert((tree->tip_count + tree->inner_count * 3) == preOrderNodes.size());
    unsigned int unique = 0;
    unsigned int duplicate = 0;
    for (auto node: preOrderNodes) {
      auto newNode = _RIDToSubtree[getRepeatIndex(node)];
      if (newNode == node) {
        unique++;
      } else {
        duplicate++;
      }
    }
    //Logger::info << "unique: " << unique << std::endl;
    //Logger::info << "dup: " << duplicate << std::endl;
  }

  void enable() {
    _enabled = true;
  }

  void disable() {
    resetCache();
    _enabled = false;
  }

  /*
   *  Get the repeat index of subtree. If this subtree is found for
   *  the first time, a new repeat index will be generated
   */
  unsigned int getRepeatIndex(pll_unode_t *subtree) {
    assert(_geneToSpecies.size());
    auto queried = _subtreeToRID.find(subtree);
    unsigned int res = -1;
    if (queried != _subtreeToRID.end()) {
      res = queried->second;
    } else {
      res = addNewSubtree(subtree);
    }
    _NIDToRID[subtree->node_index] = res;
    return res;
  }
  
  unsigned int getRepeatIndexNoCheck(pll_unode_t *subtree) {
    //std::cerr << (subtree->next ? "inner" : "leaf") << std::endl;
    assert(_NIDToRID.size() > subtree->node_index);
    auto res = _NIDToRID[subtree->node_index];
    if (res == static_cast<unsigned int>(-1)) {
      assert(!subtree->next);
      res = addNewSubtree(subtree);
      _NIDToRID[subtree->node_index] = res;
    }
    return res;
  }

  unsigned int geneToSpecies(unsigned int gene) {return _geneToSpecies[gene];}
private:
  unsigned int addNewSubtree(pll_unode_t *subtree) {
    auto newIndex = _RIDToSubtree.size();
    _RIDToSubtree.push_back(subtree);
    _subtreeToRID[subtree] = newIndex;
    return newIndex;
  }

  std::unordered_map<pll_unode_t *, unsigned int, hashing_func_subtree, key_equal_fn_subtree> _subtreeToRID;
  std::vector<pll_unode_t *> _RIDToSubtree;  // Repeat Id to subtree
  std::vector<unsigned int> _NIDToRID;  // Node Id to subtree
  std::vector<unsigned int> _geneToSpecies;
  bool _enabled;
};

