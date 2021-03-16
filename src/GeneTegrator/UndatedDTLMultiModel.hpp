#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include "ConditionalClades.hpp"
#include <maths/ScaledValue.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDTLMultiModel: public BaseReconciliationModel {
public: 
  UndatedDTLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile);

  virtual ~UndatedDTLMultiModel() {}

  virtual void setRates(const RatesVector &) {};
  virtual double computeLogLikelihood();
  virtual bool inferMLScenario(Scenario &, bool stochastic = false) {
    (void)stochastic;
    return false;}
  
private:
  
  
  ConditionalClades _ccp;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PT; // Transfer probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  double _transferExtinctionSum;

  /**
   *  All intermediate results needed to compute the reconciliation likelihood
   *  each gene node has one DTLCLV object
   *  Each DTLCLV gene  object is a function of the DTLCLVs of the direct children genes
   */
  struct DTLCLV {
    DTLCLV():
      _survivingTransferSums(REAL())
    {}

    DTLCLV(unsigned int speciesNumber):
      _uq(speciesNumber, REAL()),
      _survivingTransferSums(REAL())
    {
    }
    // probability of a gene node rooted at a species node
    std::vector<REAL> _uq;

    // sum of transfer probabilities. Can be computed only once
    // for all species, to reduce computation complexity
    REAL _survivingTransferSums;
  };
  
  std::vector<DTLCLV> _dtlclvs;



  REAL getLikelihoodFactor() const;
  virtual void recomputeSpeciesProbabilities();
  void computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba);

  void mapGenesToSpecies();
  
};


template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &geneTreesFile):
  BaseReconciliationModel(speciesTree,
      geneSpeciesMapping,
      info),
  _ccp(geneTreesFile),
  _PD(speciesTree.getNodesNumber(), 0.2),
  _PL(speciesTree.getNodesNumber(), 0.2),
  _PT(speciesTree.getNodesNumber(), 0.2),
  _PS(speciesTree.getNodesNumber(), 1.0),
  _uE(speciesTree.getNodesNumber(), 0.0)
{
  std::vector<REAL> zeros(speciesTree.getNodesNumber());
  DTLCLV nullCLV(this->_allSpeciesNodesCount);
  _dtlclvs = std::vector<DTLCLV>(2 * (_ccp.getCladesNumber()), nullCLV);
  for (unsigned int e = 0; e < _speciesTree.getNodesNumber(); ++e) {
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] /= sum;
  }
  mapGenesToSpecies();
}

template <class REAL>
double UndatedDTLMultiModel<REAL>::computeLogLikelihood()
{ 
  if (_ccp.skip()) {
    return 0.0;
  }
  beforeComputeLogLikelihood();
  onSpeciesTreeChange(nullptr);
  for (CID cid = 0; cid < _ccp.getCladesNumber(); ++cid) {
    auto &clv = _dtlclvs[cid];
    auto &uq = clv._uq;
    auto sum = REAL();
    for (auto speciesNode: _allSpeciesNodes) {
      computeProbability(cid, 
          speciesNode, 
          uq[speciesNode->node_index]);
      sum += uq[speciesNode->node_index];
    }
    sum /= this->_allSpeciesNodes.size();
    clv._survivingTransferSums = sum;
  }
  auto rootCID = _ccp.getCladesNumber() - 1;
  REAL res = REAL();
  for (auto speciesNode: _allSpeciesNodes) {
    res += _dtlclvs[rootCID]._uq[speciesNode->node_index];
  }
  // the root correction makes sure that UndatedDTLMultiModel and
  // UndatedDTL model are equivalent when there is one tree per
  // family: the UndatedDTLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDTL, so we multiply back here
  REAL rootCorrection(double(_ccp.getRootsNumber()));
  return log(res) - log(getLikelihoodFactor()) + log(rootCorrection);
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::recomputeSpeciesProbabilities()
{
  for (auto speciesNode: _allSpeciesNodes) {
    _uE[speciesNode->node_index] = REAL(0.0);
  }
  _transferExtinctionSum = REAL();
  for (unsigned int it = 0; it < 4; ++it) {
    for (auto speciesNode: _allSpeciesNodes) {
      auto e = speciesNode->node_index;
      double proba(_PL[e]);
      proba += _uE[e] * _uE[e] * _PD[e];
      proba += _transferExtinctionSum * _PT[e] * _uE[e];
      if (this->getSpeciesLeft(speciesNode)) {
        auto leftSpid = this->getSpeciesLeft(speciesNode)->node_index;
        auto rightSpid = this->getSpeciesRight(speciesNode)->node_index;
        proba += _uE[leftSpid]  * _uE[rightSpid] * _PS[e];
      }
      _uE[e] = proba;
    }
    _transferExtinctionSum = 0.0;
    for (auto speciesNode: _allSpeciesNodes) {
      _transferExtinctionSum += _uE[speciesNode->node_index];
    }
    _transferExtinctionSum /= this->_allSpeciesNodes.size();
  }
}


template <class REAL>
void UndatedDTLMultiModel<REAL>::computeProbability(CID cid, 
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
 
  // for internal gene nodes
  for (const auto &cladeSplit: _ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left; 
    auto cidRight = cladeSplit.right;
    REAL splitProba = REAL();
    if (not isSpeciesLeaf) {
      // S event;
      temp = _dtlclvs[cidLeft]._uq[f] * _dtlclvs[cidRight]._uq[g] * _PS[e]; 
      scale(temp);
      splitProba += temp;
      temp = _dtlclvs[cidRight]._uq[f] * _dtlclvs[cidLeft]._uq[g] * _PS[e]; 
      scale(temp);
      splitProba += temp;
    }
    // D event
    temp = _dtlclvs[cidLeft]._uq[e] * _dtlclvs[cidRight]._uq[e] * _PD[e];
    scale(temp);
    splitProba += temp;

    // T event
    
    temp =  _dtlclvs[cidLeft]._survivingTransferSums * _PT[e];
    temp *= _dtlclvs[cidRight]._uq[e];
    scale(temp);
    splitProba += temp;
   
    temp =  _dtlclvs[cidRight]._survivingTransferSums * _PT[e];
    temp *= _dtlclvs[cidLeft]._uq[e];
    scale(temp);
    splitProba += temp;
   
    // add this split contribution to the total
    proba += splitProba * cladeSplit.frequency;
  }
  if (not isSpeciesLeaf) {
    // SL event
    temp = _dtlclvs[cid]._uq[f] * (_uE[g] * _PS[e]);
    scale(temp);
    proba += temp;
    temp = _dtlclvs[cid]._uq[g] * (_uE[f] * _PS[e]);
    scale(temp);
    proba += temp;
  }
}
  
template <class REAL>
REAL UndatedDTLMultiModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}


  template <class REAL>
void UndatedDTLMultiModel<REAL>::mapGenesToSpecies()
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


