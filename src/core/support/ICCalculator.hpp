#pragma once
#include <vector>
#include <string>
#include <IO/Families.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <util/types.hpp>

// counts[spid][famid][gid] 
using IntersectionCounts = 
std::vector<std::vector<std::vector<unsigned int> > >;

using QuadriCounts = 
std::vector<std::vector<UInt3> >;


class ICCalculator {
public:
  ICCalculator(const std::string &referenceTreePath,
      const Families &families,
      bool paralogy);
  
  
  
private:
  // trees
  PLLRootedTree _rootedReferenceTree;
  PLLUnrootedTree _referenceTree;
  pll_unode_t *_referenceRoot;
  PerCoreGeneTrees _perCoreGeneTrees;
  TaxaSet _allSPID;
  std::vector<unsigned int> _speciesSubtreeSizes;
  unsigned int _taxaNumber;
  std::vector<unsigned int> _refNodeIndexToBranchIndex;
  unsigned int _maxBranchIndex;
 
  // parameters
  bool _paralogy;

  // intermediate results
  IntersectionCounts _interCounts;
  QuadriCounts _quadriCounts;

  // evaluation trees data
  std::vector<std::unique_ptr<PLLUnrootedTree> > _evaluationTrees;
 
  // debug
  std::vector<std::string> _spidToString;

  std::vector<double> _qpic; 
  std::vector<double> _eqpic; 
  std::vector<double> _localPP; 

  // _isPartitionInFamily[famid][spid]
  std::vector<std::vector<bool> > _isPartitionInFamily;

  void _readTrees();
  void _computeIntersections();
  void _computeQuadriCounts();

  unsigned int _getQuadripartitionCount(
    unsigned int famid,
    const std::array<unsigned int, 4> &refQuadriparition,
    const std::array<unsigned int, 3> &evalTripartition);
  unsigned int _getQuadripartitionCountPro(
    unsigned int famid,
    const std::array<unsigned int, 4> &refQuadriparition,
    const std::array<unsigned int, 3> &evalTripartition);
  void _computeRefBranchIndices();
  std::string _getNewickWithScore(std::vector<double> &branchScores, const std::string &scoreName);

};
