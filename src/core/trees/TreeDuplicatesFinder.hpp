#pragma once

#include <unordered_map>
#include <vector>

typedef struct pll_unode_s pll_unode_t;
typedef std::unordered_map<pll_unode_t*, unsigned int> TreeDuplicates;

class PerCoreGeneTrees;

class TreeDuplicatesFinder {
public:
  static void findDuplicates(PerCoreGeneTrees &perCoreGeneTrees, TreeDuplicates &duplicates);
};


template <typename Value> 
class SubtreeCache
{
public:
  SubtreeCache(const TreeDuplicates &duplicates);
  
  void resetAll();
  void resetNode(); 

  bool isValueComputed(pll_unode_t *key) const {return _isValid[getIndex(key)];}
  
  void setValue(pll_unode_t *key, const Value &value) {
    auto index = getIndex(key);
    _values[index] = value;
    _isValid[index] = true;
  }
  const Value &getValue(pll_unode_t *key) const { 
    auto index = getIndex(key);
    assert(_isValid[index]);
    return _values[index];
  } 

private:
  unsigned int getIndex(pll_unode_t *key) const {return _duplicates.at(key);}
  const TreeDuplicates &_duplicates;
  std::vector<Value> _values;
  std::vector<bool> _isValid;
};


