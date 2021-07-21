#ifdef CODE_DISABLED
#pragma once

#include <ccp/SpeciesSplits.hpp>
#include <trees/PLLRootedTree.hpp>

class SpeciesSplits;
class PLLRootedTree;

using BID = unsigned int; 
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


class SpeciesSplitScore {
public:
  SpeciesSplitScore(PLLRootedTree &speciesTree,
      const SpeciesSplits &splits);
  

private:
  std::vector<unsigned int> _bidToNodeIndex;
  std::vector<BID> _nodeIndexToBid;
  BranchSet _emptyBranchSet;
  // _path[leaf->node_inde
  std::vector<std::vector<BranchSet> > _paths;
  std::unordered_map<std::string, unsigned int> _labelToSpid;
  
  
  void fillPathsRec(unsigned int fromSpid, 
      pll_rnode_t *node,
      pll_rnode_t* previousNode,
      BranchSet &path);
};
#endif
