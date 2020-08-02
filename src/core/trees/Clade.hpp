#pragma once

#include <set>
#include <cstddef>
#include <iostream>
#include <unordered_map>
#include <util/enums.hpp>


class PLLUnrootedTree;
class PLLRootedTree;
class GeneSpeciesMapping;
using CladeSet = std::set<size_t>;

class Clade {
public:
  Clade();

  void mergeWith(const Clade &clade);
  void addId(unsigned int id);
  size_t getHash() const {return _hash;}
  friend std::ostream& operator<<(std::ostream& os, const Clade &c) 
  {
    os << "[(";
    for (auto id: c._ids) {
      os << id << ",";
    } 
    os << ")" << c._hash << "]";
    return os;
  }
  
  Clade getComplement(Clade &allTaxa);

  static Clade getMaximumClade(PLLRootedTree &tree);

  static CladeSet buildCladeSet(PLLRootedTree &tree);

  static CladeSet buildCladeSet(PLLUnrootedTree &tree,
      const GeneSpeciesMapping &geneSpeciesMapping,
      const StringToUintMap &speciesLabelToInt);
private:
  void _recomputeHash();
  
  std::set<unsigned int> _ids;
  size_t _hash;
};

struct Counter
{
  struct value_type { template<typename T> value_type(const T&) { } };
  void push_back(const value_type&) { ++count; }
  size_t count = 0;
};

  template<typename T1, typename T2>
size_t intersection_size(const T1& s1, const T2& s2)
{
  Counter c;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(c));
  return c.count;
}

