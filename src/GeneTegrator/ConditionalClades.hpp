#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include "bitvector.hpp"

using CID = unsigned int;

using CCPClade = genesis::utils::Bitvector;
using CladeCounts = std::unordered_map<CID, unsigned int>;
using SubcladeCounts = std::vector<CladeCounts>;
using OrderedClades = std::set<CCPClade>;

struct CladeSplit {
  CID parent;
  CID left;
  CID right;
  double frequency;
  CladeSplit():
    frequency(0.0)
  {
  }
};

using CladeSplits = std::vector<CladeSplit>;
using CladeToCID = std::unordered_map<CCPClade, CID>;
using CIDToClade = std::vector<CCPClade>;
using CIDToLeaf = std::unordered_map<unsigned int, std::string>;
class ConditionalClades {
public:
  ConditionalClades(const std::string &newickFile);
  void printContent() const; 
  unsigned int getCladesNumber() const {return _cladeToCID.size();}
  const CladeSplits &getCladeSplits(CID cid) const {return _allCladeSplits[cid];} 
  bool isLeaf(CID cid) const;
  std::string getLeafLabel(CID cid) const;
  unsigned int getRootsNumber() const;
  unsigned int getInputTreesNumber() const {return _inputTrees;}
  unsigned int getUniqueInputTreesNumber() const {return _uniqueInputTrees;}
  const std::unordered_map<unsigned int, std::string> &
    getCidToLeaves() {return _CIDToLeaf;}
  bool skip() const {return _skip;} 
private:
  bool _skip;
  std::vector<std::string> _idToLeaf;
  CIDToLeaf _CIDToLeaf;
  std::vector<CladeSplits> _allCladeSplits;
  CladeToCID _cladeToCID;
  CIDToClade _CIDToClade;
  unsigned int _inputTrees;
  unsigned int _uniqueInputTrees;
private:
  void _fillCCP(CladeCounts &cladeCounts,
      SubcladeCounts &subcladeCounts);
};



