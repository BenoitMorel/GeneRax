#pragma once

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <ccp/ConditionalClades.hpp>
#include <utility>

class PLLUnrootedTree;
class GeneSpeciesMapping;

using Split = std::pair<CID, CID>; 

struct PairHash
{
  template <class T1, class T2>
  std::size_t operator() (const std::pair<T1, T2> &pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

using SplitCounts = std::unordered_map<Split, unsigned int, PairHash>;

class SpeciesSplits {
public:
  SpeciesSplits(const std::unordered_set<std::string> &speciesLabels);

  void addGeneTree(PLLUnrootedTree &geneTree,
      const GeneSpeciesMapping &mapping);
  
  void addGeneTree(const std::string &newickFile,
      const GeneSpeciesMapping &mapping);

  unsigned int distinctSplitsNumber() const {return _splitCounts.size();}
  unsigned int nonDistinctSplitsNumber() const;

  const SplitCounts &getSplitCounts() const {return _splitCounts;}

private: 
  std::vector<std::string> _spidToLabel;
  std::unordered_map<std::string, unsigned int> _labelToSpid;
  unsigned int _speciesNumber;
  CladeToCID _cladeToCid;
  CIDToClade _cidToClade;
  std::vector<unsigned int> _cidToSpeciesNumber;
  SplitCounts _splitCounts;

  void _treatNodeSplit(unsigned int nodeId,
      unsigned int leftNodeId,
      unsigned int rightNodeId,
      std::vector<CID> &geneNodeToCid);
  void _addSplit(Split &split);
};


