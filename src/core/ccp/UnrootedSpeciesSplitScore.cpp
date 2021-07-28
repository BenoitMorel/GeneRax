#include "UnrootedSpeciesSplitScore.hpp"
#include <ccp/UnrootedSpeciesSplitScore.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>

static const unsigned int INVALID_BID = static_cast<unsigned int>(-1);
void UnrootedSpeciesSplitScore::fillPathsRec(unsigned int fromSpid, 
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

UnrootedSpeciesSplitScore::UnrootedSpeciesSplitScore(PLLRootedTree &rootedSpeciesTree,
    const SpeciesSplits &splits):
  _splits(splits),
  _speciesTree(std::make_unique<PLLUnrootedTree>(rootedSpeciesTree)),
  _bidToNodeIndex(_speciesTree->getLeavesNumber() - 3),
  _nodeIndexToBid(_speciesTree->getDirectedNodesNumber(), INVALID_BID),
  _emptyBranchSet(_speciesTree->getLeavesNumber() - 3),
  _labelToSpid(splits.getLabelToSpid())
{
  std::vector<BranchSet> emptySets(_speciesTree->getLeavesNumber(), _emptyBranchSet);
  _paths.resize(_speciesTree->getLeavesNumber());
  std::fill(_paths.begin(), _paths.end(), emptySets);
}
  
void UnrootedSpeciesSplitScore::updateSpeciesTree(PLLRootedTree &rootedSpeciesTree)
{
  _speciesTree = std::make_unique<PLLUnrootedTree>(rootedSpeciesTree);
  // map internal branches to branch ids
  unsigned int spid = 0;
  for (auto node: _speciesTree->getBranchesDeterministic()) {
    if (node->next && node->back->next) {
      _bidToNodeIndex[spid] = (node->node_index);
      _nodeIndexToBid[node->node_index] = spid;
      _nodeIndexToBid[node->back->node_index] = spid;
      spid++;
    }
  }
  
  // compute branch paths between each pair of species leaves
  for (auto leaf: _speciesTree->getLeaves()) {
    BranchSet path = _emptyBranchSet;
    fillPathsRec(_labelToSpid[leaf->label], leaf->back->next->back, path);
    fillPathsRec(_labelToSpid[leaf->label], leaf->back->next->next->back, path);
  }

}

double UnrootedSpeciesSplitScore::getScoreHack()
{
  static const bool weightCladeSize = false;//true; // does not matter that much
  static const bool weightBranchNumber = true; // important
  static const bool useDenominator = true;
  static const double branchPow = 2.0;
  // compute the score
  std::vector<double> score(_bidToNodeIndex.size(), 0.0);
  std::vector<double> denominator(_bidToNodeIndex.size(), 0.0);
  const auto &counts = _splits.getSplitCounts();
  auto begin = ParallelContext::getBegin(counts.size());
  auto end = ParallelContext::getEnd(counts.size());
  for (unsigned int i = begin; i < end; ++i) {
    auto cid1 = counts[i].first.first;
    auto cid2 = counts[i].first.second;
    const auto &clade1 = _splits.getClade(cid1);
    const auto &clade2 = _splits.getClade(cid2);
    auto count = counts[i].second;
    auto leftBranchSet = getRelevantBranches(cid1);
    auto rightBranchSet = getRelevantBranches(cid2);
    auto branchIntersection = leftBranchSet & rightBranchSet;
    double splitWeight = count;
    if (weightCladeSize) {
      splitWeight *= (clade1.count() + clade2.count());
    }
    if (branchIntersection.count() == 0) {
      // disjoint branch sets: the clade agrees with the species tree
      auto link = getAnyBranchPath(_splits.getClade(cid1), _splits.getClade(cid2));
      auto branchUnion = leftBranchSet | rightBranchSet;
      auto nonUnion = ~branchUnion;
      auto compatibleBranches = link & nonUnion;
      if (weightBranchNumber) {
        splitWeight /= pow(double(compatibleBranches.count()), branchPow);
      }
      for (unsigned int bid = 0; bid < compatibleBranches.size(); ++bid) {
        if (compatibleBranches[bid]) {
          score[bid] += splitWeight;
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
    ParallelContext::sumDouble(score[i]);
    ParallelContext::sumDouble(denominator[i]);
    auto ratio = score[i] / denominator[i];
    res += log(0.0001 + ratio);
  }
  return res;
}

double UnrootedSpeciesSplitScore::getScore()
{
  return getScoreHack();
  static const bool weightCladeSize = false;//true; // does not matter that much
  static const bool weightBranchNumber = true; // important
  static const bool useDenominator = true;
  static const double branchPow = 2.0;
  // compute the score
  std::vector<double> score(_bidToNodeIndex.size(), 0.0);
  std::vector<double> denominator(_bidToNodeIndex.size(), 0.0);
  const auto &counts = _splits.getSplitCounts();
  auto begin = ParallelContext::getBegin(counts.size());
  auto end = ParallelContext::getEnd(counts.size());
  for (unsigned int i = begin; i < end; ++i) {
    auto cid1 = counts[i].first.first;
    auto cid2 = counts[i].first.second;
    const auto &clade1 = _splits.getClade(cid1);
    const auto &clade2 = _splits.getClade(cid2);
    auto count = counts[i].second;
    auto leftBranchSet = getRelevantBranches(cid1);
    auto rightBranchSet = getRelevantBranches(cid2);
    auto branchIntersection = leftBranchSet & rightBranchSet;
    double splitWeight = count;
    if (weightCladeSize) {
      splitWeight *= (clade1.count() + clade2.count());
    }
    if (branchIntersection.count() == 0) {
      // disjoint branch sets: the clade agrees with the species tree
      auto link = getAnyBranchPath(_splits.getClade(cid1), _splits.getClade(cid2));
      auto branchUnion = leftBranchSet | rightBranchSet;
      auto nonUnion = ~branchUnion;
      auto compatibleBranches = link & nonUnion;
      if (weightBranchNumber) {
        splitWeight /= pow(double(compatibleBranches.count()), branchPow);
      }
      for (unsigned int bid = 0; bid < compatibleBranches.size(); ++bid) {
        if (compatibleBranches[bid]) {
          score[bid] += splitWeight;
          denominator[bid] += splitWeight;
        }
      } 
    } else {
      // non empty intersection: the clade disagrees with the species tree
      for (unsigned int bid = 0; bid < branchIntersection.size(); ++bid) {
        if (branchIntersection[bid]) {
          denominator[bid] += splitWeight;
        }
      }
    }
  }
  double res = 0;
  for (unsigned int i = 0; i < _bidToNodeIndex.size(); ++i) {
    ParallelContext::sumDouble(score[i]);
    ParallelContext::sumDouble(denominator[i]);
    if (useDenominator) {
      if (denominator[i] != 0.0) {
        auto ratio = score[i] / denominator[i];
        res += log(0.0001 + ratio);
      }
    } else  {
      auto ratio = score[i];
      res += log(0.0001 + ratio);
    }
  }
  return res;
}


BranchSet UnrootedSpeciesSplitScore::getRelevantBranches(unsigned int cid)
{
  // todo: this could be precomputed in a much more efficient way
  // (traverse the clades recursively to build the branch sets)
  BranchSet res = _emptyBranchSet;
  const auto &clade = _splits.getClade(cid);
  unsigned int previousBid = INVALID_BID;
  for (unsigned int bid = 0; bid < clade.size(); ++bid) {
    if (clade[bid]) {
      if (previousBid != INVALID_BID) {
        res |= _paths[previousBid][bid];
      }
      previousBid = bid;
    }
  }
  return res; 
}

BranchSet UnrootedSpeciesSplitScore::getAnyBranchPath(const CCPClade &c1, const CCPClade &c2)
{
  unsigned int bid1 = 0;
  unsigned int bid2 = 0;
  while (!c1[bid1]) { bid1++; }
  while (!c2[bid2]) { bid2++; }
  return _paths[bid1][bid2];
}


