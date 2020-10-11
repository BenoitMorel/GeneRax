#pragma once
#include <vector>
#include <string>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <util/types.hpp>




class ICCalculatorSlow {
public:
  ICCalculatorSlow(const std::string &referenceTreePath,
      const Families &families,
      bool paralogy);

  void printNQuartets(unsigned int n);
private:
  // reference tree data
  PLLRootedTree _rootedReferenceTree;
  PLLUnrootedTree _referenceTree;
  TaxaSet _allSPID;
  unsigned int _taxaNumber;
  std::vector<unsigned int> _refNodeIndexToBranchIndex;
  unsigned int _maxBranchIndex;

  // evaluation trees data
  std::vector<std::unique_ptr<PLLUnrootedTree> > _evaluationTrees;
  
  // quartet information
  std::vector<SPID> _quartetCounts;


  // scores
  // scores indexed with branch indices
  std::vector<double> _lqic;
  std::vector<double> _qpic; 
  std::vector<double> _eqpic; 


  // debug
  std::vector<std::string> _spidToString;

  void _readTrees(const Families &families);
  void _computeQuartets(bool paralogy);
  void _computeQuartetsForTree(PLLUnrootedTree &tree);
  void _computeQuartetsForTreePro(PLLUnrootedTree &tree);
  void _initScores();
  void _computeScores();
  void _processNodePair(pll_unode_t *u, pll_unode_t *v);


  void _computeRefBranchIndices();
  unsigned int _getLookupIndex(unsigned int a, 
      unsigned int b,
      unsigned int c,
      unsigned int d) {
    return a + _taxaNumber * (b + _taxaNumber * (c + _taxaNumber * d));
  }
  
  double _getQic(SPID a, SPID b, SPID c, SPID d);
  void _getSpidUnderNode(pll_unode_t *node, TaxaSet &taxa);

  void _printQuartet(SPID a, SPID b, SPID c, SPID d);
  void _updateFromMetaquartet(const MetaQuartet &m);
  std::string _getNewickWithScore(std::vector<double> &branchScores, const std::string &scoreName);
};


