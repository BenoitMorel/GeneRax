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
  /**
   *  @param referenceTreePath 
   *  @param families information about evaluation trees
   *  @param eqpicRadius max size of the branch path to consider for eqpic computation
   *  @param paralogy if set to true, only consider speciation-driven quartets
   */
  ICCalculator(const std::string &referenceTreePath,
      const Families &families,
      int eqpicRadius,
      bool paralogy);

  /*
   *  Export computed scores in newick format
   *  SupportTriplets corresponds to the three alternative
   *  quartet topology frequencies.
   *  Support corresponds to the frequency of the quartets
   *  present in the species tree
   *
   *  @param outputQPIC filepath to the qpic newick output
   *  @param outputEQPIC filepath to the eqpic newick output
   *  @param outputSupport filepath to the support newick output
   *  @param outputSupportTriplets filepath to the qpic newick output
   *
   */
  void exportScores(const std::string &outputQPIC,
      const std::string &outputEQPIC,
      const std::string &outputSupport,
      const std::string &outputSupportTriplets);

  static void computeScores(PLLRootedTree &tree,
      const Families &families,
      bool paralogy,
      int eqpicRadius,
      const std::string &tempPrefix,
      std::vector<double> &idToSupport);

private:
  // trees
  PLLRootedTree _rootedReferenceTree;
  PLLUnrootedTree _referenceTree;
  corax_unode_t *_referenceRoot;
  PerCoreGeneTrees _perCoreGeneTrees;
  TaxaSet _allSPID;
  std::vector<unsigned int> _speciesSubtreeSizes;
  unsigned int _taxaNumber;
  std::vector<unsigned int> _refNodeIndexToBranchIndex;
  unsigned int _maxBranchIndex;
 
  // parameters
  bool _paralogy;
  int _eqpicRadius;
  // intermediate results
  IntersectionCounts _interCounts;
  QuadriCounts _quadriCounts;

  // evaluation trees data
  std::vector<std::unique_ptr<PLLUnrootedTree> > _evaluationTrees;
 
  // debug
  std::vector<std::string> _spidToString;

  std::vector<double> _qpic; 
  std::vector<double> _eqpic; 
  std::vector<double> _localSupport1; 
  std::vector<double> _localSupport2; 
  std::vector<double> _localSupport3; 
  
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
  std::string _getNewickWithThreeScores();
  std::string _getNewickWithScore(std::vector<double> &branchScores, const std::string &scoreName);


};
