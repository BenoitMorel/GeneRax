#pragma once

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <ccp/ConditionalClades.hpp>
#include <utility>
#include <map>

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

using SplitCountsMap = std::map<Split, unsigned int>;
//using SplitCountsMap = std::unordered_map<Split, unsigned int, PairHash>;
using SplitCount = std::pair<Split, unsigned int>;
using SplitCountVector = std::vector<SplitCount>;

class SpeciesSplits {
public:
  SpeciesSplits(const std::unordered_set<std::string> &speciesLabels,
      bool acceptTrivialClade);

  void addGeneTree(PLLUnrootedTree &geneTree,
      const GeneSpeciesMapping &mapping);
  
  void addGeneTree(const std::string &newickFile,
      const GeneSpeciesMapping &mapping);

  void computeVector();

  unsigned int distinctSplitsNumber() const {return _splitCountVector.size();}
  unsigned int nonDistinctSplitsNumber() const;

  const SplitCountVector &getSplitCounts() const {return _splitCountVector;}
  const CIDToClade &getClades() const {return _cidToClade;}
  const CCPClade &getClade(CID cid) const {return _cidToClade[cid];}
  const std::unordered_map<std::string, unsigned int> &getLabelToSpid() const {return _labelToSpid;}
private:
  bool _acceptTrivialClade;
  std::vector<std::string> _spidToLabel;
  std::unordered_map<std::string, unsigned int> _labelToSpid;
  unsigned int _speciesNumber;
  CladeToCID _cladeToCid;
  CIDToClade _cidToClade;
  std::vector<unsigned int> _cidToSpeciesNumber;
  SplitCountsMap _splitCountsMap;
  SplitCountVector _splitCountVector;

  void _treatNodeSplit(unsigned int nodeId,
      unsigned int leftNodeId,
      unsigned int rightNodeId,
      std::vector<CID> &geneNodeToCid);
  void _addSplit(Split &split);
};


