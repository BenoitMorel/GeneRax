#include "ICCalculator.hpp"

#include <array>
#include <algorithm>
#include <unordered_map>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/LibpllParsers.hpp>
#include <trees/DSTagger.hpp>

ICCalculator::ICCalculator(const std::string &referenceTreePath,
      const Families &families,
      bool paralogy):
  _rootedReferenceTree(referenceTreePath),
  _referenceTree(_rootedReferenceTree),
  _taxaNumber(0)
{
  _readTrees(families);
  _computeRefBranchIndices();
  _computeIntersections();
  _computeQuadriCounts();
  Logger::info << _getNewickWithScore(_qpic, std::string("QPIC")) << std::endl;
  Logger::info << _getNewickWithScore(_eqpic, std::string("EPIC")) << std::endl;
}

void ICCalculator::_readTrees(const Families &families)
{
  Logger::timed << "[IC computation] Reading trees..." << std::endl; 
  std::unordered_map<std::string, SPID> speciesLabelToSpid;
  
  for (auto speciesLabel: _referenceTree.getLabels()) {
    auto spid = static_cast<SPID>(speciesLabelToSpid.size());
    speciesLabelToSpid[speciesLabel] = spid;
    _allSPID.insert(spid);
    _spidToString.push_back(speciesLabel);
  }
  for (auto &leaf: _referenceTree.getLeaves()) {
    leaf->clv_index = speciesLabelToSpid[std::string(leaf->label)]; 
  }
  _taxaNumber = _allSPID.size();
  for (auto &family: families) {
    GeneSpeciesMapping mappings;
    mappings.fill(family.mappingFile, family.startingGeneTree);
    auto evaluationTree = std::make_unique<PLLUnrootedTree>(family.startingGeneTree);
    for (auto &leaf: evaluationTree->getLeaves()) {
      leaf->clv_index = speciesLabelToSpid[mappings.getSpecies(std::string(leaf->label))];
    }
    _evaluationTrees.push_back(std::move(evaluationTree));
  }
}

void fillWithChildren(pll_unode_t *node, TaxaSet &set)
{
  if (!node->next) {
    set.insert(node->clv_index);
  } else {
    fillWithChildren(node->next->back, set);
    fillWithChildren(node->next->next->back, set);
  }
}

static unsigned int intersectionSize(const TaxaSet &s1, const TaxaSet &s2)
{
  if (s1.size() > s2.size()) {
    return intersectionSize(s2, s1);
  }
  unsigned int res = 0;
  for (auto element: s1) {
    if (s2.find(element) != s2.end()) {
      res++;
    }
  }
  return res;
}

void ICCalculator::_computeIntersections()
{
  Logger::timed << "[IC Computation] Precomputing intersections..." <<std::endl;
  auto speciesNodeCount = _referenceTree.getDirectedNodesNumber();
  auto familyCount= _evaluationTrees.size();
  _interCounts.clear();
  _interCounts.resize(familyCount);
  std::vector<unsigned int> speciesZeros(speciesNodeCount, 0);
  for (unsigned int famid = 0; famid < _evaluationTrees.size(); ++famid) {
    _interCounts[famid].resize(
        _evaluationTrees[famid]->getDirectedNodesNumber());
    for (auto &geneCounts:_interCounts[famid]) {
      geneCounts = speciesZeros;
    }
  }
  std::vector<TaxaSet> speciesSets(speciesNodeCount);
  for (auto speciesNode: _referenceTree.getPostOrderNodes()) {
    auto spid = speciesNode->node_index;
    fillWithChildren(speciesNode, speciesSets[spid]);
  }
  for (unsigned int famid = 0; famid < _evaluationTrees.size(); ++famid) {
    if (famid % 100 == 0) {
      Logger::timed << famid << "/" << _evaluationTrees.size() << std::endl;
    }
    auto &geneTree = _evaluationTrees[famid];
    for (auto geneNode: geneTree->getPostOrderNodes()) {
      TaxaSet geneSet;
      auto geneid = geneNode->node_index;
      // TODO: this could be more efficient (see Astral III paper)
      fillWithChildren(geneNode, geneSet);
      for (unsigned int spid = 0; spid < speciesNodeCount; ++spid) { 
        // TODO: this could also be more efficient
        auto interSize = intersectionSize(speciesSets[spid], geneSet);
        _interCounts[famid][geneid][spid] = interSize;
      }
    }
  }
}

unsigned int ICCalculator::_getQuadripartitionCount(
    unsigned int famid,
    const std::array<unsigned int, 4> &refQuadriparition,
    const std::array<unsigned int, 3> &evalTripartition)
{
  auto &T = evalTripartition;
  auto &M = refQuadriparition;
  const auto &q = _interCounts[famid];
  unsigned int res = 0;
  auto A = M[0];
  auto B = M[1];
  auto C = M[2];
  auto D = M[3];
  for (unsigned int i = 0; i < 3; ++i) { // three tripartitions
    auto AiBi = q[T[i]][A] * q[T[i]][B];
    auto CiDi = q[T[i]][C] * q[T[i]][D];
    for (unsigned int j = 0; j < 3; ++j) {
      if (i == j) {
        continue;
      }
      unsigned int k = (3 - (i + j));
      assert(k != i && k != j);
      auto CjDk = q[T[j]][C] * q[T[k]][D];
      auto AjBk = q[T[j]][A] * q[T[k]][B];
      res += AiBi * CjDk;
      res += CiDi * AjBk;
    }
  }
  return res;
}
   

/**
 *  u and v must be oriented toward each other
 *  and must be internal nodes
 *  Returns one of the three possible quadripartitions
 *  defined by u  and v. 
 *  quartetIndex == 0 -> "real" topology AB|CD
 *  quartetIndex == 1 -> alternative topology AC|BD
 *  quartetIndex == 2 -> alternative topology AD|BC
 */
static UInt4 getQuadripartition(pll_unode_t *u,
    pll_unode_t *v,
    unsigned int quartetIndex)
{
  auto A = u->next->back->node_index;
  auto B = u->next->next->back->node_index;
  auto C = v->next->back->node_index;
  auto D = v->next->next->back->node_index;
  UInt4 res;
  switch (quartetIndex) {
  case 0:
    res[0] = A;
    res[1] = B;
    res[2] = C;
    res[3] = D;
    break;
  case 1:
    res[0] = A;
    res[1] = C;
    res[2] = B;
    res[3] = D;
    break;
  case 2:
    res[0] = A;
    res[1] = D;
    res[2] = B;
    res[3] = C;
    break;
  default:
    assert(false);
  }
  return res;
}
    
static UInt3 getTripartition(pll_unode_t *u)
{
  auto A = u->back->node_index;
  auto B = u->next->back->node_index;
  auto C = u->next->next->back->node_index;
  return UInt3({A, B, C});
}

static double getLogScore(const std::array<unsigned int, 3> &q) 
{
  if (q[0] == 0 && q[1] == 0 && q[2] == 0) {
    return 0.0;
  }
  auto qsum = std::accumulate(q.begin(), q.end(), 0);
  double qic = 1.0;
  for (unsigned int i = 0; i < 3; ++i) {
    if (q[i] != 0) {
      auto p = double(q[i]) / double(qsum);
		  qic += p * log(p) / log(3);
    }
  }
  if (q[0] >= q[1] && q[0] >= q[1]) {
    return qic;
  } else {
    return -qic;
  }
}


void ICCalculator::_computeQuadriCounts()
{
  Logger::timed << "Computing scores..." << std::endl;
  auto branchNumbers = _referenceTree.getLeavesNumber() * 2 - 3;
  _qpic = std::vector<double>(branchNumbers, 1.0);
  _eqpic = std::vector<double>(branchNumbers, 1.0);
  unsigned int speciesNodeCount = _referenceTree.getDirectedNodesNumber();
  auto familyCount = _evaluationTrees.size();
  _quadriCounts.clear();
  _quadriCounts.resize(speciesNodeCount);
  for (auto &ucount: _quadriCounts) {
    ucount.resize(speciesNodeCount);
    for (auto &vcount: ucount) {
      vcount = {0, 0, 0};
    }
  }
  std::vector<pll_unode_t *> speciesInnerNodes;
  for (auto node: _referenceTree.getInnerNodes()) {
    speciesInnerNodes.push_back(node);
  }
  
  std::vector<std::vector<UInt3> > tripartitions(familyCount);
  for (unsigned int famid = 0; famid < familyCount; ++famid) {
    auto &geneTree = _evaluationTrees[famid];
    for (auto geneNode: geneTree->getPostOrderNodes(true)) {
      tripartitions[famid].push_back(getTripartition(geneNode));
    }
  }
  


  for (unsigned int i = 0; i < speciesInnerNodes.size(); ++i) {
    Logger::timed << i <<"/" << speciesInnerNodes.size() << std::endl;
    auto unode = speciesInnerNodes[i];
    for (unsigned int j = i+1; j < speciesInnerNodes.size(); ++j) {
      std::array<unsigned int, 3> counts = {0, 0, 0};
      auto vnode = speciesInnerNodes[j];
      assert(unode->next && vnode->next);
      assert(unode->next != vnode && unode->next->next != vnode);
      assert(vnode->next != unode && vnode->next->next != unode);
      std::vector<pll_unode_t *> branchPath;
      _referenceTree.orientTowardEachOther(&unode, 
          &vnode,
          branchPath);
      std::vector<unsigned int> branchIndices;
      for (auto branch: branchPath) {
        branchIndices.push_back(_refNodeIndexToBranchIndex[branch->node_index]);
      }
      std::array<UInt4, 3> quadripartitions;
      for (unsigned int topology = 0; topology < 3; ++topology) {
        quadripartitions[topology] = getQuadripartition(unode, vnode, topology);
      }
      for (unsigned int famid = 0; famid < familyCount; ++famid) {
        for (auto &tripartition: tripartitions[famid]) {
          for (unsigned int topology = 0; topology < 3; ++topology) {
            counts[topology] += _getQuadripartitionCount(famid,
                quadripartitions[topology], 
                tripartition);
          }
        }
      }
      auto qpic = getLogScore(counts);
      if (branchIndices.size() == 1) {
        auto branchIndex = branchIndices[0];
        _qpic[branchIndex] = qpic;
      }
      for (auto branchIndex: branchIndices) {
        _eqpic[branchIndex] = std::min(_eqpic[branchIndex], qpic);
      }
    }
  }
}

std::string ICCalculator::_getNewickWithScore(std::vector<double> &branchScores, const std::string &scoreName)
{
  for (auto node: _referenceTree.getPostOrderNodes()) {
    if (!node->next || !node->back->next) {
      continue;
    }
    auto branchIndex = _refNodeIndexToBranchIndex[node->node_index];
    auto score = branchScores[branchIndex];
    auto scoreStr = scoreName + std::string(" = ") + std::to_string(score);
    free(node->label);
    node->label = (char *)calloc(scoreStr.size() + 1, sizeof(char));
    strcpy(node->label, scoreStr.c_str());
  }
  return _referenceTree.getNewickString(); 
}


void ICCalculator::_computeRefBranchIndices()
{
  Logger::timed << "[IC computation] Assigning branch indices..." << std::endl; 
  unsigned int currBranchIndex = 0;
  auto branches = _referenceTree.getBranches();
  _refNodeIndexToBranchIndex.resize(branches.size() * 2);
  std::fill(_refNodeIndexToBranchIndex.begin(),
      _refNodeIndexToBranchIndex.end(),
      static_cast<unsigned int>(-1));
  for (auto &branch: branches) {
    _refNodeIndexToBranchIndex[branch->node_index] = currBranchIndex;
    _refNodeIndexToBranchIndex[branch->back->node_index] = currBranchIndex;
    currBranchIndex++;
  }
  assert(currBranchIndex == _taxaNumber * 2 - 3);
  for (auto v: _refNodeIndexToBranchIndex) {
    assert(v != static_cast<unsigned int>(-1));
  }
}

