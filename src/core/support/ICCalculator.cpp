#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include "ICCalculator.hpp"

#include <array>
#include <algorithm>
#include <unordered_map>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/LibpllParsers.hpp>
#include <trees/DSTagger.hpp>
#include <sstream>
#include <cmath>
#include <maths/incbeta.h>
#include <IO/ParallelOfstream.hpp>

ICCalculator::ICCalculator(const std::string &referenceTreePath,
      const Families &families,
      bool paralogy):
  _rootedReferenceTree(referenceTreePath),
  _referenceTree(_rootedReferenceTree),
  _referenceRoot(_referenceTree.getVirtualRoot(_rootedReferenceTree)),
  _perCoreGeneTrees(families),
  _taxaNumber(0),
  _paralogy(paralogy)
{
}

void ICCalculator::computeScores(const std::string &outputQPIC,
    const std::string &outputEQPIC,
    const std::string &outputSupport)
{
  _readTrees();
  _computeRefBranchIndices();
  _computeIntersections();
  _computeQuadriCounts();

  Logger::info << "Writing species tree with QPIC scores in " << outputQPIC << std::endl;
  Logger::info << "Writing species tree with EQPIC scores in " << outputEQPIC << std::endl;
  Logger::info << "Writing species tree with quartet support scores in " << outputSupport << std::endl;
  ParallelOfstream osQPIC(outputQPIC);
  osQPIC << _getNewickWithScore(_qpic, std::string()) << std::endl;
  ParallelOfstream osEQPIC(outputEQPIC);
  osEQPIC << _getNewickWithScore(_eqpic, std::string()) << std::endl;
  ParallelOfstream osSupport(outputSupport);
  osSupport << _getNewickWithScore(_localPP, std::string()) << std::endl;
}


void ICCalculator::_readTrees()
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
  for (auto &geneTreeDesc: _perCoreGeneTrees.getTrees()) {
    auto &mappings = geneTreeDesc.mapping;
    auto evaluationTree = geneTreeDesc.geneTree;
    for (auto &leaf: evaluationTree->getLeaves()) {
      leaf->clv_index = speciesLabelToSpid[mappings.getSpecies(std::string(leaf->label))];
    }
    _evaluationTrees.push_back(
        std::unique_ptr<PLLUnrootedTree>(evaluationTree));
    geneTreeDesc.ownTree = false;
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
  _speciesSubtreeSizes.resize(_referenceTree.getDirectedNodesNumber());
  for (auto speciesNode: _referenceTree.getPostOrderNodes()) {
    auto spid = speciesNode->node_index;
    fillWithChildren(speciesNode, speciesSets[spid]);
    _speciesSubtreeSizes[spid] = speciesSets[spid].size();
  }
  for (unsigned int famid = 0; famid < _evaluationTrees.size(); ++famid) {
    auto &geneTree = _evaluationTrees[famid];
    
    std::unique_ptr<DSTagger> tagger(nullptr);
    if (_paralogy) {
      tagger = std::make_unique<DSTagger>(*geneTree);
    }
    for (auto geneNode: geneTree->getPostOrderNodes()) {
      auto geneid = geneNode->node_index;
      if (_paralogy && tagger->isDuplication(geneid))  {
        continue;
      }
      // TODO: this could be more efficient (see Astral III paper)
      TaxaSet geneSet;
      if (_paralogy && tagger->goesUp(geneNode)) {
        tagger->fillUpTraversal(geneNode, geneSet);
      } else {
        fillWithChildren(geneNode->back, geneSet);
      }
      for (unsigned int spid = 0; spid < speciesNodeCount; ++spid) { 
        auto interSize = intersectionSize(speciesSets[spid], geneSet);
        _interCounts[famid][geneid][spid] = interSize;
      }
    }
  }
}

unsigned int ICCalculator::_getQuadripartitionCountPro(
    unsigned int famid,
    const std::array<unsigned int, 4> &refQuadriparition,
    const std::array<unsigned int, 3> &evalTripartition)
{
  unsigned int res = 0;
  auto A = 0;
  auto B = 1;
  auto C = 2;
  auto D = 3;
  /*
  auto &T = evalTripartition;
  auto &M = refQuadriparition;
  const auto &q = _interCounts[famid];
  if (!q[T[0]][A] && !q[T[1]][A] && !q[T[2]][A]) {
    return 0;
  }
  */
  unsigned int v[3][4];
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 4; ++j) {
      v[i][j] = _interCounts[famid][evalTripartition[i]][refQuadriparition[j]];
    }
  }
  auto x1a = v[0][A]; 
  auto x1b = v[0][B]; 
  auto x1c = v[0][C]; 
  auto x1d = v[0][D]; 
  auto x2a = v[1][A]; 
  auto x2b = v[1][B]; 
  auto x2c = v[1][C]; 
  auto x2d = v[1][D]; 
  auto ya = v[2][A]; 
  auto yb = v[2][B]; 
  auto yc = v[2][C]; 
  auto yd = v[2][D]; 

  // balanced anchors
  res += x1a*x1b*x2c*x2d + x1c*x1d*x2a*x2b + x1a*x1b*(x2c*yd+yc*x2d) + x1c*x1d*(x2a*yb+ya*x2b) + x2a*x2b*(x1c*yd+yc*x1d) + x2c*x2d*(x1a*yb+ya*x1b);
  return res;
}
   
unsigned int ICCalculator::_getQuadripartitionCount(
    unsigned int famid,
    const std::array<unsigned int, 4> &refQuadriparition,
    const std::array<unsigned int, 3> &evalTripartition)
{
  unsigned int res = 0;
  auto A = 0;
  auto B = 1;
  auto C = 2;
  auto D = 3;
  /*
  auto &T = evalTripartition;
  auto &M = refQuadriparition;
  const auto &q = _interCounts[famid];
  if (!q[T[0]][A] && !q[T[1]][A] && !q[T[2]][A]) {
    return 0;
  }
  */
  unsigned int v[3][4];
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 4; ++j) {
      v[i][j] = _interCounts[famid][evalTripartition[i]][refQuadriparition[j]];
    }
  }
  for (unsigned int i = 0; i < 3; ++i) { // three tripartitions
    auto AiBi = v[i][A] * v[i][B];
    auto CiDi = v[i][C] * v[i][D];
    auto j = (i + 1) % 3;
    auto k = (i + 2) % 3;
    auto CjDk = v[j][C] * v[k][D];
    auto CkDj = v[k][C] * v[j][D];
    auto AjBk = v[j][A] * v[k][B];
    auto AkBj = v[k][A] * v[j][B];
    res += AiBi * (CjDk + CkDj);
    res += CiDi * (AjBk + AkBj);
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
    
UInt3 getTripartition(pll_unode_t *u,
    DSTagger *tagger)
{
  if (tagger) {
    tagger->orientUp(u);
  }
  auto y = u->node_index;
  auto x1 = u->next->node_index;
  auto x2 = u->next->next->node_index;
  return UInt3({x1, x2, y});
}

static double getLogScore(const std::array<unsigned long, 3> &q) 
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
static double computeH(double z, double zsum, double lambda)
{
  auto term1 = std::beta(z + 1.0, zsum - z  + 2.0 * lambda);
  auto term2 = incbeta(z + 1.0, zsum - z + 2.0 * lambda, 1.0/3.0);
  return term1 * (1.0 - term2);
}

static double computeLocalPP(const std::array<unsigned long, 3> &counts,
    unsigned long ratio)
{
  // We use notations from astral-pro paper
  std::array<double, 3> z;
  for (unsigned int i = 0; i < 3; ++i) {
    z[i] = double(counts[i]) / double(ratio);
  }

  std::array<double, 3> h;
  double lambda = 0.5; // value from astral-pro paper
  auto zsum = std::accumulate(z.begin(), z.end(), 0);
  for (unsigned int i = 0; i < 3; ++i) {
    h[i] = computeH(z[i], zsum, lambda);
  }
  return z[0] / zsum;
  double den = h[0];
  den += pow(2.0, z[1] - z[0]) * h[1];
  den += pow(2.0, z[2] - z[0]) * h[2];
  return h[0] / den;
}


void ICCalculator::_computeQuadriCounts()
{
  Logger::timed << "Computing scores..." << std::endl;
  auto branchNumbers = _referenceTree.getLeavesNumber() * 2 - 3;
  _qpic = std::vector<double>(branchNumbers, 1.0);
  _eqpic = std::vector<double>(branchNumbers, 1.0);
  _localPP = std::vector<double>(branchNumbers, 1.0);
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
    std::unique_ptr<DSTagger> tagger(nullptr);
    if (_paralogy) {
      tagger = std::make_unique<DSTagger>(*geneTree);
    }
    for (auto geneNode: geneTree->getInnerNodes()) { 
      tripartitions[famid].push_back(getTripartition(geneNode, tagger.get())); 
    }
  }


  for (unsigned int i = 0; i < speciesInnerNodes.size(); ++i) {
    auto unode = speciesInnerNodes[i];
    for (unsigned int j = 0; j < speciesInnerNodes.size(); ++j) {
    //for (unsigned int j = i+1; j < speciesInnerNodes.size(); ++j) {
      if (i == j) {
        continue;
      }
      std::array<unsigned long, 3> counts = {0, 0, 0};
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
            if(_paralogy) {
              counts[topology] += _getQuadripartitionCountPro(famid,
                  quadripartitions[topology], 
                  tripartition);
            } else {
              counts[topology] += _getQuadripartitionCount(famid,
                  quadripartitions[topology], 
                  tripartition);
            }
          }
        }
      }
      if (branchIndices.size() == 1) {
        for (unsigned int topology = 0; topology < 3; ++topology) {
          ParallelContext::sumULong(counts[topology]);
        }
      }
      auto qpic = getLogScore(counts);
      if (branchIndices.size() == 1) {
        auto branchIndex = branchIndices[0];
        _qpic[branchIndex] = qpic;
        unsigned long ratio = 1;
        ratio *= _speciesSubtreeSizes[unode->next->back->node_index];
        ratio *= _speciesSubtreeSizes[unode->next->next->back->node_index];
        ratio *= _speciesSubtreeSizes[vnode->next->back->node_index];
        ratio *= _speciesSubtreeSizes[vnode->next->next->back->node_index];
        _localPP[branchIndex] = computeLocalPP(counts, ratio);
      }
      for (auto branchIndex: branchIndices) {
        _eqpic[branchIndex] = std::min(_eqpic[branchIndex], qpic);
      }
    }
  }
}


struct ScorePrinter
{
  ScorePrinter(const std::string& prefix,
      const std::vector<double> &score,
      const std::vector<unsigned int> &refNodeIndexToBranchIndex):
    prefix(prefix),
    score(score),
    refNodeIndexToBranchIndex(refNodeIndexToBranchIndex)
  {}
  std::string prefix;
  const std::vector<double> &score;
  const std::vector<unsigned int> &refNodeIndexToBranchIndex;
  void operator()(pll_unode_t *node, 
      std::stringstream &ss)
  {
    if (node->next && node->back->next) {
      auto idx = refNodeIndexToBranchIndex[node->node_index];
      ss << prefix << score[idx];
    } else {
      if (!node->next && node->label) {
        ss << node->label;
      }
    }
    ss << ":" << node->length;
  }
};

std::string ICCalculator::_getNewickWithScore(std::vector<double> &branchScores, const std::string &scoreName)
{
  ScorePrinter printer(scoreName.size() ? (scoreName  + " = ") : "", 
      branchScores, 
      _refNodeIndexToBranchIndex);
  return _referenceTree.getNewickString(printer,
      _referenceRoot,
      true); 
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

