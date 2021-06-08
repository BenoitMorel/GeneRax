#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include "bitvector.hpp"
#include <iostream>

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
  ConditionalClades() {}
  ConditionalClades(const std::string &inputFile, 
      bool fromBinary = false);
  void printContent() const; 
  void printStats() const;
  unsigned int getCladesNumber() const {return _cladeToCID.size();}
  const CladeSplits &getCladeSplits(CID cid) const {return _allCladeSplits[cid];} 
  bool isLeaf(CID cid) const;
  std::string getLeafLabel(CID cid) const;
  unsigned int getRootsNumber() const;
  unsigned int getInputTreesNumber() const {return _inputTrees;}
  unsigned int getUniqueInputTreesNumber() const {return _uniqueInputTrees;}
  const std::unordered_map<unsigned int, std::string> &
    getCidToLeaves() {return _CIDToLeaf;}
  bool skip() const {return false;}
  //bool skip() const {return  _uniqueInputTrees == _inputTrees;} 

  void serialize(const std::string &outputFile);
  void unserialize(const std::string &inputFile);
private:
  unsigned int _inputTrees;
  unsigned int _uniqueInputTrees;
  std::vector<std::string> _idToLeaf;
  CIDToLeaf _CIDToLeaf;
  CladeToCID _cladeToCID;
  CIDToClade _CIDToClade;
  std::vector<CladeSplits> _allCladeSplits;
private:
  void _fillCCP(CladeCounts &cladeCounts,
      SubcladeCounts &subcladeCounts);
};



