#include "ICCalculator.hpp"

#include <array>
#include <unordered_map>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>


ICCalculator::ICCalculator(const std::string &referenceTreePath,
      const Families &families):
  _rootedReferenceTree(referenceTreePath),
  _referenceTree(_rootedReferenceTree),
  _evaluationTrees()
{
  _readTrees(referenceTreePath, families);
 
  /*
  for (auto u: _referenceTree.getInnerNodes()) {
    for (auto v: _referenceTree.getInnerNodes()) {
      _processNodePair(u, v);
    }
  }
  */
}


void ICCalculator::_readTrees(const std::string &referenceTreePath,
    const Families &families)
{
  std::unordered_map<std::string, SPID> labelToId;
  for (auto p: _rootedReferenceTree.getLabelToIntMap()) {
    labelToId[p.first] = static_cast<SPID>(p.second);
    _allTaxa.insert(labelToId[p.first]);
  }
  for (auto &leaf: _rootedReferenceTree.getLeaves()) {
    leaf->data = (void*)labelToId[std::string(leaf->label)];
  }
  for (auto &leaf: _referenceTree.getLeaves()) {
    leaf->data = (void*)labelToId[std::string(leaf->label)];
  }

  for (auto &family: families) {
    GeneSpeciesMapping mappings;
    mappings.fill(family.mappingFile, family.startingGeneTree);
    auto evaluationTree = std::make_unique<PLLUnrootedTree>(family.startingGeneTree);
    for (auto &leaf: evaluationTree->getLeaves()) {
      leaf->data = (void*)labelToId[mappings.getSpecies(std::string(leaf->label))];
    }
    _evaluationTrees.push_back(std::move(evaluationTree));
  }

}
  
void ICCalculator::_processNodePair(pll_rnode_t *u, pll_rnode_t *v)
{
  assert(u->left && u->right);
  assert(v->left && v->right);
   
}

