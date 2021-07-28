#pragma once

#include <ccp/SpeciesSplits.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <ccp/SpeciesSplitScore.hpp>

class SpeciesSplits;
class PLLRootedTree;

using BranchSet = genesis::utils::Bitvector;


/**
 *  The BID is an unique identifer for an internal branch in an unrooted tree. 
 *  If the tree is rooted, the root and its right branches are not mapped to
 *  a BID. 
 *
 *  A path between two leaves is the (unordered) set of internal 
 *  branches (terminal branches are excluded) from one leaf to another. 
 *  A path is represented with a BranchSet (a set of BIDs)
 *
 *
 */


class UnrootedSpeciesSplitScore: public SpeciesSplitScore {
public:
  UnrootedSpeciesSplitScore(PLLRootedTree &rootedSpeciesTree,
      const SpeciesSplits &splits);
  virtual ~UnrootedSpeciesSplitScore() {}  
  virtual void updateSpeciesTree(PLLRootedTree &rootedSpeciesTree);
  virtual double getScore();
private:
  double getScoreHack();
  const SpeciesSplits &_splits;
  std::unique_ptr<PLLUnrootedTree> _speciesTree;
  std::vector<unsigned int> _bidToNodeIndex;
  std::vector<BID> _nodeIndexToBid;
  BranchSet _emptyBranchSet;
  std::vector<std::vector<BranchSet> > _paths;
  std::unordered_map<std::string, unsigned int> _labelToSpid;
  void fillPathsRec(unsigned int fromSpid, 
    pll_unode_t *node,
    BranchSet &path);
  
  // returns all the branchs that are located between at least
  // one pair of leaves in the clade identified by cid
  BranchSet getRelevantBranches(unsigned int cid);
  // select one leaf from c1 and one leaf from c2, and returns
  // the path between those two leaves
  BranchSet getAnyBranchPath(const CCPClade &c1, const CCPClade &c2);
};
