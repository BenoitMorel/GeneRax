#pragma once

#include <unordered_map>
#include <likelihoods/LibpllEvaluation.hpp>
#include <vector>


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
  }

  void resetCache() {
    _subtreeToRID.clear();
    _RIDToSubtree.clear();
  }

  

  /*
   *  Get the repeat index of subtree. If this subtree is found for
   *  the first time, a new repeat index will be generated
   */
  unsigned int getRepeatIndex(pll_unode_t *subtree) {
    assert(_geneToSpecies.size());
    auto queried = _subtreeToRID.find(subtree);
    if (queried != _subtreeToRID.end()) {
      return queried->second;
    }
    return addNewSubtree(subtree);
  }
  
  unsigned int getRepeatIndexNoCheck(pll_unode_t *subtree) {
    assert(_geneToSpecies.size());
    auto queried = _subtreeToRID.find(subtree);
    assert(queried != _subtreeToRID.end());
    return queried->second;
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
  std::vector<pll_unode_t *> _RIDToSubtree;
  std::vector<unsigned int> _geneToSpecies;
};

