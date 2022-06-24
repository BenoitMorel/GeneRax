#include "RootedSpeciesSplitScore.hpp"
#include <ccp/RootedSpeciesSplitScore.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>

static const SPID INVALID_SPID = static_cast<unsigned int>(-1);
static const unsigned int INVALID_NODE_INDEX = static_cast<unsigned int>(-1);

RootedSpeciesSplitScore::RootedSpeciesSplitScore(PLLRootedTree &speciesTree,
    const SpeciesSplits &splits):
  _splits(splits),
  _speciesTree(&speciesTree),
  _spidToNodeIndex(_speciesTree->getNodesNumber()),
  _nodeIndexToSpid(_speciesTree->getNodesNumber(), INVALID_SPID)
  //_labelToSpid(splits.getLabelToSpid())
{
}
  
void RootedSpeciesSplitScore::updateSpeciesTree(PLLRootedTree &speciesTree)
{
  _speciesTree = &speciesTree;
  _speciesTree->buildLCACache();
  // map internal branches to branch ids
  SPID maxSpid = 0;
  for (auto node: _speciesTree->getLeaves()) {
    std::string label(node->label);
    auto spid = _splits.getLabelToSpid().at(label);
    _spidToNodeIndex[spid] = (node->node_index);
    _nodeIndexToSpid[node->node_index] = spid;
    maxSpid = std::max(spid, maxSpid);
  }
  for (auto node: _speciesTree->getInnerNodes()) {
    maxSpid++;
    _spidToNodeIndex[maxSpid] = (node->node_index);
    _nodeIndexToSpid[node->node_index] = maxSpid;
  }
}


corax_rnode_t *RootedSpeciesSplitScore::getLCA(unsigned int cid)
{
  const auto &clade = _splits.getClade(cid);
  corax_rnode_t *lca = nullptr;
  for (unsigned int spid = 0; spid < clade.size(); ++spid) {
    if (clade[spid]) {
      if (lca == nullptr) {
        lca = _speciesTree->getNode(_spidToNodeIndex[spid]);
      } else {
        lca = _speciesTree->getLCA(lca->node_index, _spidToNodeIndex[spid]);
      }
    }
  } 
  return lca;
}


double RootedSpeciesSplitScore::getScore()
{
  // compute the score
  std::vector<double> score(_spidToNodeIndex.size(), 0.0);
  std::vector<double> denominator(_spidToNodeIndex.size(), 0.0);
  const auto &counts = _splits.getSplitCounts();
  auto begin = ParallelContext::getBegin(counts.size());
  auto end = ParallelContext::getEnd(counts.size());
  for (unsigned int i = begin; i < end; ++i) {
    auto cid1 = counts[i].first.first;
    auto cid2 = counts[i].first.second;
    auto count = counts[i].second;
    auto lca1 = getLCA(cid1);
    auto lca2 = getLCA(cid2);
    auto lca = _speciesTree->getLCA(lca1, lca2);
    auto lcaSpid = _nodeIndexToSpid[lca->node_index];
    if (_speciesTree->areParents(lca1, lca2)) {
      // MISMATCH
      denominator[lcaSpid] += count; 
    } else {
      // MATCH
      score[lcaSpid] += count;
      denominator[lcaSpid] += count;
    }
  }
  double res = 0;
  double den = 1.0;
  for (unsigned int spid = 0; spid < _spidToNodeIndex.size(); ++spid) {
    auto node = _speciesTree->getNode(_spidToNodeIndex[spid]);
    if (!node->left) {
      continue;
    }
    ParallelContext::sumDouble(score[spid]);
    ParallelContext::sumDouble(denominator[spid]);
    res += score[spid];
  }
  return res / den;
}
