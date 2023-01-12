#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <maths/bitvector.hpp>
#include <iostream>
#include <util/enums.hpp>

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
  double deviation;
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
      CCPRooting ccpRooting);
  void printContent() const; 
  void printStats() const;
  unsigned int getCladesNumber() const {return _cladeToCID.size();}
  const CladeSplits &getCladeSplits(CID cid) const {return _allCladeSplits[cid];} 
  bool isLeaf(CID cid) const;
  std::string getLeafLabel(CID cid) const;
  unsigned int getRootsNumber() const;
  unsigned int getInputTreesNumber() const {return _inputTrees;}
  unsigned int getUniqueInputTreesNumber() const {return _uniqueInputTrees;}
  unsigned int getLeafNumber() const {return _CIDToLeaf.size();}
  const std::unordered_map<unsigned int, std::string> &
    getCidToLeaves() {return _CIDToLeaf;}
  bool skip() const {return false;}
  //bool skip() const {return  _uniqueInputTrees == _inputTrees;} 

  void serialize(const std::string &outputFile);
  void unserialize(const std::string &inputFile);

  bool madRooting() const {return _ccpRooting == CCPRooting::MAD;}

private:
  unsigned int _inputTrees;
  unsigned int _uniqueInputTrees;
  CCPRooting _ccpRooting;
  std::vector<std::string> _idToLeaf;
  CIDToLeaf _CIDToLeaf;
  CladeToCID _cladeToCID;
  CIDToClade _CIDToClade;
  std::vector<CladeSplits> _allCladeSplits;
private:
  void _fillCCP(CladeCounts &cladeCounts,
      SubcladeCounts &subcladeCounts,
      std::unordered_map<unsigned int, double> *CIDToDeviation = nullptr);
};



