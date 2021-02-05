#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <set>

using CCPClade = std::vector<bool>;
using CladeCounts = std::unordered_map<CCPClade, unsigned int>;
using SubcladeCounts = std::unordered_map<CCPClade, CladeCounts>;
using OrderedClades = std::set<CCPClade>;

class ConditionalClades {
public:
  ConditionalClades(const std::string &newickFile);
  
  std::set<CCPClade>::reverse_iterator begin() {
    return _internalClades.rbegin();
  }

  std::set<CCPClade>::reverse_iterator end() {
    return _internalClades.rend();
  }
  
  void printContent(); 

private:
  std::vector<std::string> _idToLeaf;

  CladeCounts _cladeCounts;
  SubcladeCounts _subcladeCounts;
  OrderedClades _internalClades;
};



