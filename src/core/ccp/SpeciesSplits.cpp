#include "SpeciesSplits.hpp"
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLUnrootedTree.hpp>

const unsigned int INVALID_NODE_ID = static_cast<unsigned int>(-1);

SpeciesSplits::SpeciesSplits(const std::unordered_set<std::string> &speciesLabels):
  _speciesNumber(speciesLabels.size())
{
  for (auto &label: speciesLabels) {
    unsigned int spid = _labelToSpid.size();
    _labelToSpid.insert({label, spid});
    _spidToLabel.push_back(label);
    CCPClade clade(_speciesNumber);
    clade.set(spid);
    _cidToClade.push_back(clade);
    _cidToSpeciesNumber.push_back(1);
    _cladeToCid.insert({clade, spid});
  }
}
  
void SpeciesSplits::_treatNodeSplit(unsigned int nodeId,
    unsigned int leftNodeId,
    unsigned int rightNodeId, 
    std::vector<CID> &geneNodeToCid)
{
  auto leftCid = geneNodeToCid[leftNodeId];
  auto rightCid = geneNodeToCid[rightNodeId];
  auto &leftClade = _cidToClade[leftCid];
  auto &rightClade = _cidToClade[rightCid];
  bool isSpeciation = (leftClade & rightClade).count() == 0;
  auto clade = leftClade | rightClade;
  CID cid;
  auto cidIt = _cladeToCid.find(clade);
  if (_cladeToCid.find(clade) == _cladeToCid.end()) {
    cid = _cladeToCid.size();
    _cladeToCid.insert({clade, cid});
    _cidToClade.push_back(clade);
    _cidToSpeciesNumber.push_back(clade.count());
  } else {
    cid = cidIt->second;
  }
  if (nodeId != INVALID_NODE_ID) {
    geneNodeToCid[nodeId] = cid;
  }
  if (isSpeciation) {
    Split split;
    split.first = leftCid;
    split.second = rightCid;
    _addSplit(split);
  }
}


void SpeciesSplits::addGeneTree(PLLUnrootedTree &geneTree,
    const GeneSpeciesMapping &mapping)
{
  auto postOrderGeneNodes = geneTree.getPostOrderNodes();
  std::vector<CID> geneNodeToCid(postOrderGeneNodes.size());
  for (auto node: postOrderGeneNodes) {
    if (!node->next) {
      // for leaves, spid == cid
      auto species = mapping.getSpecies(node->label);
      auto spid = _labelToSpid[species];
      geneNodeToCid[node->node_index] = spid;
    } else {
      auto leftNodeId = node->next->back->node_index;
      auto rightNodeId = node->next->next->back->node_index;
      _treatNodeSplit(node->node_index,
          leftNodeId, 
          rightNodeId,
          geneNodeToCid);
    }
  }
  for (auto branch: geneTree.getBranches()) {
    if (branch->next && branch->back->next) { // for all internal branches
      _treatNodeSplit(INVALID_NODE_ID,
          branch->node_index,
          branch->back->node_index,
          geneNodeToCid);
    }
  }
}

  
void SpeciesSplits::addGeneTree(const std::string &newickFile,
      const GeneSpeciesMapping &mapping)
{
  PLLUnrootedTree geneTree(newickFile);
  addGeneTree(geneTree, mapping);
}
  
unsigned int SpeciesSplits::nonDistinctSplitsNumber() const
{
  unsigned int res = 0;
  for (auto splitCount: getSplitCounts()) {
    res += splitCount.second;
  }
  return res;
}

void SpeciesSplits::_addSplit(Split &split)
{
  // Add non-trivial splits
  if (_cidToSpeciesNumber[split.first] < 2) {
    return;
  }
  if (_cidToSpeciesNumber[split.second] < 2) {
    return;
  }
  if (split.first > split.second) {
    std::swap<CID>(split.first, split.second);
  }
  auto it = _splitCounts.find(split);
  if (it != _splitCounts.end()) {
    it->second++;
  } else {
    _splitCounts.insert({split, 1});
  }
}

