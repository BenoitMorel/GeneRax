#pragma once

#include <ccp/SpeciesSplits.hpp>
#include <trees/PLLRootedTree.hpp>
#include <ccp/SpeciesSplitScore.hpp>



class RootedSpeciesSplitScore: public SpeciesSplitScore {
public:
  RootedSpeciesSplitScore(PLLRootedTree &speciesTree,
      const SpeciesSplits &splits);
  virtual ~RootedSpeciesSplitScore() {}
  virtual void updateSpeciesTree(PLLRootedTree &speciesTree);
  virtual double getScore();
private:
  const SpeciesSplits &_splits;
  PLLRootedTree *_speciesTree;
  std::vector<unsigned int> _spidToNodeIndex;
  std::vector<SPID> _nodeIndexToSpid;
  corax_rnode_t *getLCA(unsigned int cid);

};
