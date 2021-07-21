#include "SpeciesSplitScore.hpp"
#include <ccp/SpeciesSplitScore.hpp>

static const unsigned int INVALID_BID = static_cast<unsigned int>(-1);
void SpeciesSplitScore::fillPathsRec(unsigned int fromSpid, 
    pll_unode_t *node,
    BranchSet &path)
{
  if (!node->next) { 
    // we reached the end of a path
    unsigned int toSpid = _labelToSpid[node->label];
    _paths[fromSpid][toSpid] = path;
    return;
  }
  auto bid = _nodeIndexToBid[node->node_index];
  path.set(bid);
  fillPathsRec(fromSpid, node->next->back, path);
  fillPathsRec(fromSpid, node->next->next->back, path);
  path.unset(bid);
}

SpeciesSplitScore::SpeciesSplitScore(PLLRootedTree &rootedSpeciesTree,
    const SpeciesSplits &splits):
  _splits(splits),
  _speciesTree(rootedSpeciesTree),
  _nodeIndexToBid(_speciesTree.getDirectedNodesNumber(), INVALID_BID),
  _emptyBranchSet(_speciesTree.getLeavesNumber() - 3),
  _labelToSpid(splits.getLabelToSpid())
{
  // map internal branches to branch ids
  for (auto node: _speciesTree.getBranches()) {
    if (node->next) {
      unsigned int spid = _bidToNodeIndex.size();
      _bidToNodeIndex.push_back(node->node_index);
      _nodeIndexToBid[node->node_index] = spid;
      _nodeIndexToBid[node->back->node_index] = spid;
    }
  }
  std::vector<BranchSet> emptySets(_speciesTree.getLeavesNumber(), _emptyBranchSet);
  _paths.resize(_speciesTree.getLeavesNumber());
  std::fill(_paths.begin(), _paths.end(), emptySets);
  for (auto leaf: _speciesTree.getLeaves()) {
    BranchSet path = _emptyBranchSet;
    fillPathsRec(_labelToSpid[leaf->label], leaf->back->next->back, path);
    fillPathsRec(_labelToSpid[leaf->label], leaf->back->next->next->back, path);
  }
}

double SpeciesSplitScore::getScore()
{
  std::vector<double> score(_bidToNodeIndex.size(), 0.0);
  std::vector<double> denominator(_bidToNodeIndex.size(), 0.0);
  for (const auto &splitCount: _splits.getSplitCounts()) {
    auto cid1 = splitCount.first.first;
    auto cid2 = splitCount.first.second;
    auto count = splitCount.second;
    auto leftBranchSet = getRelevantBranches(cid1);
    auto rightBranchSet = getRelevantBranches(cid2);
    auto branchIntersection = leftBranchSet & rightBranchSet;
    auto splitWeight = count;
    if (branchIntersection.count() == 0) {
      // disjoint branch sets: the clade agrees with the species tree
      auto link = getAnyBranchPath(_splits.getClade(cid1), _splits.getClade(cid2));
      auto branchUnion = leftBranchSet | rightBranchSet;
      auto nonUnion = ~branchUnion;
      auto compatibleBranches = link & nonUnion;
      for (unsigned int bid = 0; bid < compatibleBranches.size(); ++bid) {
        if (compatibleBranches[bid]) {
          score[bid] += splitWeight;
          denominator[bid] += splitWeight;
        }
      } 
    } else {
      // non empty intersection: the clade disagrees with the species tree
      for (unsigned int bid = 0; bid < branchIntersection.size(); ++bid) {
        denominator[bid] += splitWeight;
      }
    }
  }

  double res = 0;
  for (unsigned int i = 0; i < _bidToNodeIndex.size(); ++i) {
    if (denominator[i] != 0.0) {
      res += score[i] / denominator[i];
    }
  }
  return res;
}


BranchSet SpeciesSplitScore::getRelevantBranches(unsigned int cid)
{
  // todo: this could be precomputed in a much more efficient way
  // (traverse the clades recursively to build the branch sets)
  BranchSet res = _emptyBranchSet;
  const auto &clade = _splits.getClade(cid);
  unsigned int previousBid = INVALID_BID;
  for (unsigned int bid = 0; bid < clade.size(); ++bid) {
    if (clade[bid]) {
      if (previousBid != INVALID_BID) {
        res = (res | _paths[previousBid][bid]);
      }
      previousBid = bid;
    }
  }
  return res; 
}

BranchSet SpeciesSplitScore::getAnyBranchPath(const CCPClade &c1, const CCPClade &c2)
{
  unsigned int bid1 = 0;
  unsigned int bid2 = 0;
  while (!c1[bid1]) { bid1++; }
  while (!c2[bid2]) { bid2++; }
  return _paths[bid1][bid2];
}


