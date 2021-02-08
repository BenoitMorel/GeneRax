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
  GeneSpeciesMapping _geneSpeciesMapping;
  ConditionalClades _ccp;
  std::vector<pll_rnode_t *> _allSpeciesNodes;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  using DLCLV = std::vector<REAL>;
  std::vector<DLCLV> _dlclvs;

  REAL getLikelihoodFactor() const;
  virtual void _recomputeSpeciesProbabilities();
  void computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba);

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
  _geneSpeciesMapping(geneSpeciesMapping),
  _ccp(geneTreesFile),
  _allSpeciesNodes(speciesTree.getPostOrderNodes()),
  _PD(speciesTree.getNodesNumber(), 0.2),
  _PL(speciesTree.getNodesNumber(), 0.2),
  _PS(speciesTree.getNodesNumber(), 1.0),
  _uE(speciesTree.getNodesNumber(), 0.0)
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
}



template <class REAL>
double UndatedDLMultiModel<REAL>::computeLogLikelihood()
{
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
  return log(res) - log(getLikelihoodFactor());
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
    std::string speciesLabel = speciesNode->label;
    std::string geneLabel = _ccp.getLeafLabel(cid);
    std::string geneSpeciesLabel = _geneSpeciesMapping.getMap().at(geneLabel);
    // todo do not compare strings...
    if (geneSpeciesLabel == speciesLabel) {
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

