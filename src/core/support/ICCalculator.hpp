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



class ICCalculator {
public:
  ICCalculator(const std::string &referenceTreePath,
      const Families &families);


private:
  PLLRootedTree _rootedReferenceTree;
  PLLUnrootedTree _referenceTree;
  std::vector<std::unique_ptr<PLLUnrootedTree> > _evaluationTrees;
  TaxaSet _allTaxa;

  void _readTrees(const std::string &referenceTreePath,
      const Families &families);
  void _processNodePair(pll_rnode_t *u, pll_rnode_t *v);
};


