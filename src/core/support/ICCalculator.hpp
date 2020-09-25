#pragma once

#include <string>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/PLLUnrootedTree.hpp>

using SPID = unsigned int;
using TaxaSet = std::unordered_set<SPID>;
struct MetaQuartet {
  std::array<TaxaSet, 4> sets;
};

using SPIDSet = std::unordered_set<SPID>;

class ICCalculator {
public:
  ICCalculator(const std::string &referenceTreePath,
      const Families &families);


private:
  // reference tree data
  PLLRootedTree _rootedReferenceTree;
  PLLUnrootedTree _referenceTree;
  TaxaSet _allTaxa;
  std::vector<SPID> _refIdToSPID;

  // evaluation trees data
  std::vector<std::unique_ptr<PLLUnrootedTree> > _evaluationTrees;
  

  void _readTrees(const std::string &referenceTreePath,
      const Families &families);
  void _mainLoop();
  void _processNodePair(pll_unode_t *u, pll_unode_t *v);
};


