#include "ICCalculator.hpp"

#include <array>
#include <unordered_map>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/LibpllParsers.hpp>

ICCalculator::ICCalculator(const std::string &referenceTreePath,
      const Families &families):
  _rootedReferenceTree(referenceTreePath),
  _referenceTree(_rootedReferenceTree),
  _evaluationTrees()
{
  _readTrees(referenceTreePath, families);
  _mainLoop();  
}


void ICCalculator::_readTrees(const std::string &referenceTreePath,
    const Families &families)
{
  std::unordered_map<std::string, SPID> labelToId;
  for (auto p: _rootedReferenceTree.getLabelToIntMap()) {
    labelToId[p.first] = static_cast<SPID>(p.second);
    _allTaxa.insert(labelToId[p.first]);
  }
  _refIdToSPID.resize(_rootedReferenceTree.getLeavesNumber());
  for (auto &leaf: _rootedReferenceTree.getLeaves()) {
    _refIdToSPID[leaf->node_index] = labelToId[std::string(leaf->label)];
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
  
void ICCalculator::_mainLoop()
{
  for (auto u: _referenceTree.getInnerNodes()) {
    for (auto v: _referenceTree.getInnerNodes()) {
      _processNodePair(u, v);
    }
  }
}
  

static SPIDSet _getSPIDClade(
    pll_unode_t *u,
    const std::vector<SPID> &idToSPID)
{
  SPIDSet set;
  for (auto id: PLLUnrootedTree::getClade(u)) {
    set.insert(idToSPID[id]);
  }
  return set;
}

void ICCalculator::_processNodePair(pll_unode_t *u, pll_unode_t *v)
{
  if (u == v) {
    return;
  }
  PLLUnrootedTree::orientTowardEachOther(&u, &v);
  assert(u != v);
  assert(u->next && v->next);
  std::array<pll_unode_t*, 4> referenceSubtrees;
  referenceSubtrees[0] = u->next->back;
  referenceSubtrees[1] = u->next->next->back;
  referenceSubtrees[2] = v->next->back;
  referenceSubtrees[3] = v->next->next->back;
  std::array<SPIDSet, 4> referenceClades;
  for (unsigned int i = 0; i < 4; ++i) {
    referenceClades[i] = _getSPIDClade(referenceSubtrees[i], 
        _refIdToSPID);
  }
}

