#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <set>

using CCPClade = std::vector<bool>;
using CladeCounts = std::unordered_map<CCPClade, unsigned int>;
using SubcladeCounts = std::unordered_map<CCPClade, CladeCounts>;
using OrderedClades = std::set<CCPClade>;

using CID = unsigned int;
struct CladeSplit {
  CID parent;
  CID left;
  CID right;
  double frequency;
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
private:
  std::vector<std::string> _idToLeaf;
  std::vector<CladeSplits> _allCladeSplits;
  CladeToCID _cladeToCID;
  CIDToClade _CIDToClade;
  
private:
  void _fillCCP(CladeCounts &cladeCounts,
      SubcladeCounts &subcladeCounts,
      OrderedClades &orderedClades);
};



