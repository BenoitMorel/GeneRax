#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include "bitvector.hpp"
//#define BOOSTIMPLEM
#define GENESISIMPLEM
#ifdef GENESISIMPLEM
using CCPClade = genesis::utils::Bitvector;
#else
using CCPClade = std::vector<bool>;
#endif
using CladeCounts = std::unordered_map<CCPClade, unsigned int>;
using SubcladeCounts = std::unordered_map<CCPClade, CladeCounts>;
using OrderedClades = std::set<CCPClade>;

using CID = unsigned int;
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

class ConditionalClades {
public:
  ConditionalClades(const std::string &newickFile);
  void printContent() const; 
  unsigned int getCladesNumber() const {return _cladeToCID.size();}
  const CladeSplits &getCladeSplits(CID cid) const {return _allCladeSplits[cid];} 
  bool isLeaf(CID cid) const;
  std::string getLeafLabel(CID cid) const;
private:
  std::vector<std::string> _idToLeaf;
  std::unordered_map<unsigned int, std::string> _CIDToLeaf;
  std::vector<CladeSplits> _allCladeSplits;
  CladeToCID _cladeToCID;
  CIDToClade _CIDToClade;
  
private:
  void _fillCCP(CladeCounts &cladeCounts,
      SubcladeCounts &subcladeCounts,
      OrderedClades &orderedClades);
};



