#pragma once

#include <unordered_map>
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

class SubtreeRepeatsCache {
public:
  SubtreeRepeatsCache(): 
    _subtreeToRID(10, hashing_func_subtree(*this), key_equal_fn_subtree(*this))
  {}

  void setGenesToSpecies(const std::vector<unsigned int> &geneToSpecies) {
    _geneToSpecies = geneToSpecies;
    resetCache();
  }

  void resetCache() {
    _subtreeToRID.clear();
    _RIDToSubtree.clear();
    _NIDToRID = std::vector<unsigned int>(_geneToSpecies.size() + 1, static_cast<unsigned int>(-1));;
  }

  pll_unode_t *getRepeat(pll_unode_t *subtree) {
    /*
    Logger::info << "getRepeat " << subtree->node_index << " ";
    if (subtree->next) {
      Logger::info << subtree->next->back->node_index << " " << subtree->next->next->back->node_index;
    }
    Logger::info << std::endl;
    */
    auto repeatIndex = getRepeatIndex(subtree);
    return _RIDToSubtree[repeatIndex];
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
};

