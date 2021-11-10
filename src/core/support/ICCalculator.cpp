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
      int eqpicRadius,
      bool paralogy):
  _rootedReferenceTree(referenceTreePath),
  _referenceTree(_rootedReferenceTree),
  _referenceRoot(_referenceTree.getVirtualRoot(_rootedReferenceTree)),
  _perCoreGeneTrees(families),
  _taxaNumber(0),
  _paralogy(paralogy),
  _eqpicRadius(eqpicRadius)
{
}

void ICCalculator::exportScores(const std::string &outputQPIC,
    const std::string &outputEQPIC,
    const std::string &outputSupport,
    const std::string &outputSupportTriplets)
{
  Logger::info << std::endl;
  Logger::timed << "[Species tree support] Read trees" << std::endl;
  _readTrees();
  _computeRefBranchIndices();
  ParallelContext::barrier();
  Logger::timed << "[Species tree support] Compute intersections" << std::endl;
  _computeIntersections();
  ParallelContext::barrier();
  Logger::timed << "[Species tree support] Compute quadricounts" << std::endl;
  _computeQuadriCounts();
  ParallelContext::barrier();
  Logger::timed << "[Species tree support] End of support computation" << std::endl;
  ParallelOfstream osQPIC(outputQPIC);
  osQPIC << _getNewickWithScore(_qpic, std::string()) << std::endl;
  ParallelOfstream osEQPIC(outputEQPIC);
  osEQPIC << _getNewickWithScore(_eqpic, std::string()) << std::endl;
  ParallelOfstream osSupport(outputSupport);
  osSupport << _getNewickWithScore(_localSupport1, std::string()) << std::endl;
  ParallelOfstream osSupportTriplet(outputSupportTriplets);
  osSupportTriplet << _getNewickWithThreeScores() << std::endl;
}


void ICCalculator::_readTrees()
{
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

void fillWithChildren(corax_unode_t *node, TaxaSet &set)
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
static UInt4 getQuadripartition(corax_unode_t *u,
    corax_unode_t *v,
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
    
UInt3 getTripartition(corax_unode_t *u,
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

void ICCalculator::_computeQuadriCounts()
{
  auto branchNumbers = _referenceTree.getLeavesNumber() * 2 - 3;
  _qpic = std::vector<double>(branchNumbers, 1.0);
  _eqpic = std::vector<double>(branchNumbers, 1.0);
  _localSupport1 = std::vector<double>(branchNumbers, 1.0);
  _localSupport2 = std::vector<double>(branchNumbers, 1.0);
  _localSupport3 = std::vector<double>(branchNumbers, 1.0);
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
  std::vector<corax_unode_t *> speciesInnerNodes;
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
    //for (unsigned int j = 0; j < speciesInnerNodes.size(); ++j) {
    for (unsigned int j = i+1; j < speciesInnerNodes.size(); ++j) {
      if (i == j) {
        continue;
      }
      std::array<unsigned long, 3> counts = {0, 0, 0};
      auto vnode = speciesInnerNodes[j];
      assert(unode->next && vnode->next);
      assert(unode->next != vnode && unode->next->next != vnode);
      assert(vnode->next != unode && vnode->next->next != unode);
      std::vector<corax_unode_t *> branchPath;
      _referenceTree.orientTowardEachOther(&unode, 
          &vnode,
          branchPath);
      if (static_cast<int>(branchPath.size()) > _eqpicRadius) {
        continue;
      }
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
      for (unsigned int topology = 0; topology < 3; ++topology) {
        ParallelContext::sumULong(counts[topology]);
      }
      auto qpic = getLogScore(counts);
      if (branchIndices.size() == 1) {
        auto branchIndex = branchIndices[0];
        _qpic[branchIndex] = qpic;
        auto sum = double(std::accumulate(counts.begin(), counts.end(), 0));
        _localSupport1[branchIndex] = double(counts[0]) / sum;
        _localSupport2[branchIndex] = double(counts[1]) / sum;
        _localSupport3[branchIndex] = double(counts[2]) / sum;
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
  void operator()(corax_unode_t *node, 
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

struct ThreeScoresPrinter
{
  ThreeScoresPrinter(const std::string& prefix,
      const std::vector<double> &score1,
      const std::vector<double> &score2,
      const std::vector<double> &score3,
      const std::vector<unsigned int> &refNodeIndexToBranchIndex):
    prefix(prefix),
    score1(score1),
    score2(score2),
    score3(score3),
    refNodeIndexToBranchIndex(refNodeIndexToBranchIndex)
  {}
  std::string prefix;
  const std::vector<double> &score1;
  const std::vector<double> &score2;
  const std::vector<double> &score3;
  const std::vector<unsigned int> &refNodeIndexToBranchIndex;
  void operator()(corax_unode_t *node, 
      std::stringstream &ss)
  {
    if (node->next && node->back->next) {
      auto idx = refNodeIndexToBranchIndex[node->node_index];
      ss << prefix
        << score1[idx] << "|"
        << score2[idx] << "| "
        << score3[idx];
    } else {
      if (!node->next && node->label) {
        ss << node->label;
      }
    }
    ss << ":" << node->length;
  }
};

std::string ICCalculator::_getNewickWithThreeScores()
{
  ThreeScoresPrinter printer("", 
      _localSupport1,
      _localSupport2,
      _localSupport3,
      _refNodeIndexToBranchIndex);
  return _referenceTree.getNewickString(printer, _referenceRoot, true); 
}

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
  
void ICCalculator::computeScores(PLLRootedTree &tree,
      const Families &families,
      bool paralogy,
      int eqpicRadius,
      const std::string &tempPrefix,
      std::vector<double> &idToSupport)
{
  bool master = ParallelContext::getRank() == 0;
  std::string tempInputTree = tempPrefix + "input.newick";
  std::string tempOutputQPIC = tempPrefix + "qpic.newick";
  std::string tempOutputEQPIC = tempPrefix + "eqpic.newick";
  std::string tempOutputSupport = tempPrefix + "support.newick";
  std::string tempOutputSupportTriplets = tempPrefix + "supportTriplets.newick";
  if  (master){
    tree.save(tempInputTree);
  }
  ParallelContext::barrier();
  ICCalculator calculator(tempInputTree, families, eqpicRadius, paralogy);
  calculator.exportScores(tempOutputQPIC, 
      tempOutputEQPIC, 
      tempOutputSupport,
      tempOutputSupportTriplets);
  ParallelContext::barrier();
  idToSupport.resize(tree.getNodesNumber());
  PLLRootedTree treeWithSupport(tempOutputEQPIC, true);
  Logger::timed << "Tree with support:" << std::endl;
  Logger::info << treeWithSupport.getNewickString() << std::endl;
  PLLRootedTree treeWithSupport2(tempOutputQPIC, true);
  Logger::timed << "Tree with support (qpic):" << std::endl;
  Logger::info << treeWithSupport2.getNewickString() << std::endl;
  auto mapping = treeWithSupport.getNodeIndexMapping(tree);
  for (auto supportedNode: treeWithSupport.getInnerNodes()) {
    auto nodeIndex = mapping[supportedNode->node_index];
    idToSupport[nodeIndex] = atof(supportedNode->label);
  }
}
    
