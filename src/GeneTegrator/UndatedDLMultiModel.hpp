#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include "ConditionalClades.hpp"
#include <maths/ScaledValue.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDLMultiModel: public BaseReconciliationModel {
public: 
  UndatedDLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile);

  virtual ~UndatedDLMultiModel() {}

  virtual void setRates(const RatesVector &rates) {};
  virtual double computeLogLikelihood();
  virtual bool inferMLScenario(Scenario &scenario, bool stochastic = false) {return false;}
  
private:
  
  
  ConditionalClades _ccp;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  using DLCLV = std::vector<REAL>;
  std::vector<DLCLV> _dlclvs;

  REAL getLikelihoodFactor() const;
  virtual void recomputeSpeciesProbabilities();
  void computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba);

  void mapGenesToSpecies();
  
};

static double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &geneTreesFile):
  BaseReconciliationModel(speciesTree,
      geneSpeciesMapping,
      info),
  _ccp(geneTreesFile),
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
  mapGenesToSpecies();
}

template <class REAL>
double UndatedDLMultiModel<REAL>::computeLogLikelihood()
{ 
  if (_ccp.skip()) {
    return 0.0;
  }
  beforeComputeLogLikelihood();
  onSpeciesTreeChange(nullptr);
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
void UndatedDLMultiModel<REAL>::recomputeSpeciesProbabilities()
{
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (getSpeciesLeft(speciesNode)) {
      c += _PS[e] * _uE[getSpeciesLeft(speciesNode)->node_index]  * 
        _uE[getSpeciesRight(speciesNode)->node_index];
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
  bool isSpeciesLeaf = !getSpeciesLeft(speciesNode);
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
    f = getSpeciesLeft(speciesNode)->node_index;
    g = getSpeciesRight(speciesNode)->node_index;
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


  template <class REAL>
void UndatedDLMultiModel<REAL>::mapGenesToSpecies()
{
  const auto &cidToLeaves = _ccp.getCidToLeaves();
  _speciesNameToId.clear();
  this->_geneToSpecies.resize(_ccp.getCladesNumber());
  for (auto node: _allSpeciesNodes) {
    if (!node->left) {
      _speciesNameToId[node->label] = node->node_index;
    }
  }
  this->_speciesCoverage = std::vector<unsigned int>(
      this->_allSpeciesNodesCount, 0);
  for (auto p: cidToLeaves) {
    auto cid = p.first;
    const auto &geneName = cidToLeaves.at(cid);
    const auto &speciesName = _geneNameToSpeciesName[geneName];
    _geneToSpecies[cid] = _speciesNameToId[speciesName];
    _speciesCoverage[_geneToSpecies[cid]]++;
  }
}

