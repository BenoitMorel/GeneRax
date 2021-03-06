#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include "ConditionalClades.hpp"
#include <maths/ScaledValue.hpp>

double log(ScaledValue v);

template <class REAL>
class UndatedDLMultiModel {
public: 
  UndatedDLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const std::string &geneTreesFile);

  virtual ~UndatedDLMultiModel() {}

  virtual double computeLogLikelihood();

private:
  
  
  PLLRootedTree &_speciesTree;
  ConditionalClades _ccp;
  std::vector<pll_rnode_t *> _allSpeciesNodes;
  unsigned int _allSpeciesNodesCount;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  using DLCLV = std::vector<REAL>;
  std::vector<DLCLV> _dlclvs;
  std::map<std::string, std::string> _geneNameToSpeciesName;
  std::unordered_map<unsigned int, unsigned int> _geneToSpecies;
  std::map<std::string, unsigned int> _speciesNameToId;

  /*
  bool _pruneSpeciesTree;
  std::vector<pll_rnode_t *> _speciesLeft;
  std::vector<pll_rnode_t *> _speciesRight;
  std::vector<pll_rnode_t *> _speciesParent;
  pll_rnode_t *_prunedRoot;
  std::vector<unsigned int> _speciesCoverage;
  */

  /*
  pll_rnode_t *getSpeciesLeft(pll_rnode_t *node) {
    return _speciesLeft[node->node_index];}
  pll_rnode_t *getSpeciesRight(pll_rnode_t *node) {
    return _speciesRight[node->node_index];}
  pll_rnode_t *getSpeciesParent(pll_rnode_t *node) {
    return _speciesParent[node->node_index];}
  pll_rnode_t *getPrunedRoot() {return _prunedRoot;}
*/

  /*
  bool fillPrunedNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd = nullptr);
*/
  REAL getLikelihoodFactor() const;
  virtual void _recomputeSpeciesProbabilities();
  void computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba);

  void initSpeciesTree();
  void mapGenesToSpecies();
  void onSpeciesTreeChange();
  
};

static double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const std::string &geneTreesFile):
  _speciesTree(speciesTree),
  _ccp(geneTreesFile),
  _allSpeciesNodes(speciesTree.getPostOrderNodes()),
  _allSpeciesNodesCount(_allSpeciesNodes.size()),
  _PD(speciesTree.getNodesNumber(), 0.2),
  _PL(speciesTree.getNodesNumber(), 0.2),
  _PS(speciesTree.getNodesNumber(), 1.0),
  _uE(speciesTree.getNodesNumber(), 0.0),
  _geneNameToSpeciesName(geneSpeciesMapping.getMap())
  //_pruneSpeciesTree(false)
{
  std::vector<REAL> zeros(speciesTree.getNodesNumber());
  _dlclvs = std::vector<std::vector<REAL> >(
      _ccp.getCladesNumber(), zeros);

  for (unsigned int e = 0; e < _speciesTree.getNodesNumber(); ++e) {
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  initSpeciesTree();
}

template <class REAL>
void UndatedDLMultiModel<REAL>::onSpeciesTreeChange()
{
  _allSpeciesNodes = _speciesTree.getPostOrderNodes();
}
  /*
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    _speciesLeft[e] = speciesNode->left;
    _speciesRight[e] = speciesNode->right;
    _speciesParent[e] = speciesNode->parent;
  }
  _prunedRoot = _speciesTree.getRoot();
  if (_pruneSpeciesTree && _speciesCoverage.size()) {
    std::vector<pll_rnode_t *> pruned(_allSpeciesNodesCount, nullptr);
    for (auto speciesNode: _allSpeciesNodes) {
      auto e = speciesNode->node_index;
      if (!speciesNode->left) {
        pruned[e] = (_speciesCoverage[e] ? speciesNode : nullptr);   
      } else {
        auto left = _speciesLeft[e];
        auto right = _speciesRight[e];
        auto prunedLeft = pruned[left->node_index];
        auto prunedRight = pruned[right->node_index];
        if (prunedLeft && prunedRight) {
          _speciesLeft[e] = prunedLeft;
          _speciesRight[e] = prunedRight;
          pruned[e] = speciesNode;
          _speciesParent[prunedLeft->node_index] = speciesNode;
          _speciesParent[prunedRight->node_index] = speciesNode;
          _prunedRoot = speciesNode;
        } else if (!prunedLeft && prunedRight) {
          pruned[e] = prunedRight;
        } else if (prunedLeft && !prunedRight) {
          pruned[e] = prunedLeft;
        }
      }
    }
    _allSpeciesNodes.clear();
    fillPrunedNodesPostOrder(getPrunedRoot(), _allSpeciesNodes);
  }
}
*/

template <class REAL>
double UndatedDLMultiModel<REAL>::computeLogLikelihood()
{ 

  double a = _ccp.getInputTreesNumber();
  double b = _ccp.getUniqueInputTreesNumber();
  if (a == b) {
    return 0.0;
  }
  onSpeciesTreeChange();
  _recomputeSpeciesProbabilities();
  for (CID cid = 0; cid < _ccp.getCladesNumber(); ++cid) {
    for (auto speciesNode: _allSpeciesNodes) {
      computeProbability(cid, 
          speciesNode, 
          _dlclvs[cid][speciesNode->node_index]);
    }
  }
  auto rootCID = _ccp.getCladesNumber() - 1;
  REAL res = REAL();
  for (auto speciesNode: _allSpeciesNodes) {
    res += _dlclvs[rootCID][speciesNode->node_index];
  }
  // the root correction makes sure that UndatedDLMultiModel and
  // UndatedDL model are equivalent when there is one tree per
  // family: the UndatedDLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDL, so we multiply back here
  REAL rootCorrection(double(_ccp.getRootsNumber()));
  return log(res) - log(getLikelihoodFactor()) + log(rootCorrection);
}

template <class REAL>
void UndatedDLMultiModel<REAL>::_recomputeSpeciesProbabilities()
{
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (speciesNode->left) {
      c += _PS[e] * _uE[speciesNode->left->node_index]  * 
        _uE[speciesNode->right->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    _uE[speciesNode->node_index] = proba;
  }
}

template <class REAL>
void UndatedDLMultiModel<REAL>::computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba
    )
{
  proba = REAL();
  bool isSpeciesLeaf = !speciesNode->left;
  auto e = speciesNode->node_index;
  if (_ccp.isLeaf(cid) && isSpeciesLeaf) {
    if (_geneToSpecies[cid] == e) {
      proba = REAL(_PS[e]);
    }
    return;
  }
  REAL temp;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  
  for (const auto &cladeSplit: _ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left; 
    auto cidRight = cladeSplit.right;
    REAL splitProba = REAL();
    if (not isSpeciesLeaf) {
      // S event;
      temp = _dlclvs[cidLeft][f] * _dlclvs[cidRight][g] * _PS[e]; 
      scale(temp);
      splitProba += temp;
      temp = _dlclvs[cidRight][f] * _dlclvs[cidLeft][g] * _PS[e]; 
      scale(temp);
      splitProba += temp;
    }
    temp = _dlclvs[cidLeft][e] * _dlclvs[cidRight][e] * _PD[e];
    scale(temp);
    splitProba += temp;
    proba += splitProba * cladeSplit.frequency;
  }
  if (not isSpeciesLeaf) {
    // SL event
    temp = _dlclvs[cid][f] * (_uE[g] * _PS[e]);
    scale(temp);
    proba += temp;
    temp = _dlclvs[cid][g] * (_uE[f] * _PS[e]);
    scale(temp);
    proba += temp;
  }
  // DL event
  proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
}
  
template <class REAL>
REAL UndatedDLMultiModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}

/*
  template <class REAL>
bool UndatedDLMultiModel<REAL>::fillPrunedNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd)
{
  bool addMyself = true;
  if (nodesToAdd) {
    addMyself = (nodesToAdd->find(node) != nodesToAdd->end());
  }
  if (getSpeciesLeft(node)) {
    assert(getSpeciesRight(node));
    addMyself |= fillPrunedNodesPostOrder(getSpeciesLeft(node), nodes, nodesToAdd);
    addMyself |= fillPrunedNodesPostOrder(getSpeciesRight(node), nodes, nodesToAdd);
  }
  if (addMyself) {
    nodes.push_back(node);
  }
  return addMyself;
}
*/
  template <class REAL>
void UndatedDLMultiModel<REAL>::initSpeciesTree()
{
  _allSpeciesNodes = _speciesTree.getPostOrderNodes();
  _allSpeciesNodesCount = _allSpeciesNodes.size();
  /*
  _speciesLeft = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesRight = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesParent = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  */
  mapGenesToSpecies();
  //onSpeciesTreeChange();
}

  template <class REAL>
void UndatedDLMultiModel<REAL>::mapGenesToSpecies()
{
  const auto &cidToLeaves = _ccp.getCidToLeaves();
  _speciesNameToId.clear();
  for (auto node: _allSpeciesNodes) {
    if (!node->left) {
      _speciesNameToId[node->label] = node->node_index;
    }
  }
  for (auto p: cidToLeaves) {
    auto cid = p.first;
    const auto &geneName = cidToLeaves.at(cid);
    const auto &speciesName = _geneNameToSpeciesName[geneName];
    _geneToSpecies[cid] = _speciesNameToId[speciesName];
  }

}

