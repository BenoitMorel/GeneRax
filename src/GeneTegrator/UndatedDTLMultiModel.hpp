#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include "MultiModel.hpp"

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDTLMultiModel: public MultiModel<REAL> {
public: 
  UndatedDTLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile);

  virtual ~UndatedDTLMultiModel() {}

  virtual void setRates(const RatesVector &);
  virtual double computeLogLikelihood();
  
private:
  
  
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
  virtual void computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell<REAL> *recCell = nullptr);
  void sampleTransferEvent(unsigned int cid,
    REAL survivingTransferSum,
    pll_rnode_t *&destSpeciesNode);

  
};


template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &geneTreesFile):
  MultiModel<REAL>(speciesTree,
      geneSpeciesMapping,
      info,
      geneTreesFile),
  _PD(speciesTree.getNodesNumber(), 0.2),
  _PL(speciesTree.getNodesNumber(), 0.2),
  _PT(speciesTree.getNodesNumber(), 0.2),
  _PS(speciesTree.getNodesNumber(), 1.0),
  _uE(speciesTree.getNodesNumber(), 0.0)
{
  std::vector<REAL> zeros(speciesTree.getNodesNumber());
  DTLCLV nullCLV(this->_allSpeciesNodesCount);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_ccp.getCladesNumber()), nullCLV);
  for (unsigned int e = 0; e < this->_speciesTree.getNodesNumber(); ++e) {
    double sum = _PD[e] + _PL[e] + _PT[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] /= sum;
  }
}
 
template <class REAL>
void UndatedDTLMultiModel<REAL>::setRates(const RatesVector &rates) 
{
  auto &dupRates = rates[0];
  auto &lossRates = rates[1];
  auto &transferRates = rates[2];
  assert(this->_allSpeciesNodesCount == dupRates.size());
  assert(this->_allSpeciesNodesCount == lossRates.size());
  assert(this->_allSpeciesNodesCount == transferRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PT = transferRates;
  _PS.resize(this->_allSpeciesNodesCount);
  for (unsigned int e = 0; e < this->_allSpeciesNodesCount; ++e) {
    if (this->_info.noDup) {
      _PD[e] = 0.0;
    }
    auto sum = _PD[e] + _PL[e] + _PT[e] + 1.0;
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] = 1.0 / sum;
  } 
  recomputeSpeciesProbabilities();
}



template <class REAL>
double UndatedDTLMultiModel<REAL>::computeLogLikelihood()
{ 
  if (this->_ccp.skip()) {
    return 0.0;
  }
  this->beforeComputeLogLikelihood();
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    auto &clv = _dtlclvs[cid];
    auto &uq = clv._uq;
    clv._survivingTransferSums = REAL();
    std::fill(uq.begin(), uq.end(), REAL());
  }  
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    auto &clv = _dtlclvs[cid];
    auto &uq = clv._uq;
    auto sum = REAL();
    for (auto speciesNode: this->_allSpeciesNodes) {
      computeProbability(cid, 
          speciesNode, 
          uq[speciesNode->node_index]);
      sum += uq[speciesNode->node_index];
    }
    sum /= this->_allSpeciesNodes.size();
    clv._survivingTransferSums = sum;
  }
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  REAL res = REAL();
  for (auto speciesNode: this->_allSpeciesNodes) {
    res += _dtlclvs[rootCID]._uq[speciesNode->node_index];
  }
  // the root correction makes sure that UndatedDTLMultiModel and
  // UndatedDTL model are equivalent when there is one tree per
  // family: the UndatedDTLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDTL, so we multiply back here
  REAL rootCorrection(double(this->_ccp.getRootsNumber()));
  auto ret = log(res) - log(getLikelihoodFactor()) + log(rootCorrection);
  return ret;
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::sampleTransferEvent(unsigned int cid,
    REAL survivingTransferSum,
    pll_rnode_t *&destSpeciesNode)
{
  REAL max = survivingTransferSum * Random::getProba();
  max *= this->_allSpeciesNodes.size();
  REAL sum = REAL();
  for (auto speciesNode: this->_allSpeciesNodes) {
    sum += _dtlclvs[cid]._uq[speciesNode->node_index];
    if (sum > max) {
      destSpeciesNode = speciesNode;
      return;
    }
  }
  assert(false);
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::recomputeSpeciesProbabilities()
{
  for (auto speciesNode: this->_allSpeciesNodes) {
    _uE[speciesNode->node_index] = REAL(0.0);
  }
  _transferExtinctionSum = REAL();
  for (unsigned int it = 0; it < 4; ++it) {
    for (auto speciesNode: this->_allSpeciesNodes) {
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
    for (auto speciesNode: this->_allSpeciesNodes) {
      _transferExtinctionSum += _uE[speciesNode->node_index];
    }
    _transferExtinctionSum /= this->_allSpeciesNodes.size();
  }
}


template <class REAL>
void UndatedDTLMultiModel<REAL>::computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell<REAL> *recCell
    )
{
  proba = REAL();
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto e = speciesNode->node_index;
  REAL maxProba = REAL();
  if (recCell) {
    recCell->event.geneNode = cid; 
    recCell->event.speciesNode = e;
    recCell->event.type = ReconciliationEventType::EVENT_None; 
    maxProba = recCell->maxProba;
  }
  if (this->_ccp.isLeaf(cid) && isSpeciesLeaf) {
    if (this->_geneToSpecies[cid] == e) {
      proba = REAL(_PS[e]);
    }
    return;
  }
  REAL temp;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
  }
 
  // for internal gene nodes
  for (const auto &cladeSplit: this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left; 
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    if (not isSpeciesLeaf) {
      // S event;
      temp = _dtlclvs[cidLeft]._uq[f] * _dtlclvs[cidRight]._uq[g] * (_PS[e] * freq); 
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        return;
      }
      temp = _dtlclvs[cidRight]._uq[f] * _dtlclvs[cidLeft]._uq[g] * (_PS[e] * freq); 
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        return;
      }
    }
    // D event
    temp = _dtlclvs[cidLeft]._uq[e] * _dtlclvs[cidRight]._uq[e] * (_PD[e] * freq);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_D;
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      return;
    }


    // T event
    
    temp =  _dtlclvs[cidLeft]._survivingTransferSums * (_PT[e] * freq);
    temp *= _dtlclvs[cidRight]._uq[e];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      sampleTransferEvent(cidLeft, 
          _dtlclvs[cidLeft]._survivingTransferSums,
          recCell->event.pllDestSpeciesNode);
      recCell->event.destSpeciesNode = 
        recCell->event.pllDestSpeciesNode->node_index;
      recCell->event.leftGeneIndex = cidRight; 
      recCell->event.rightGeneIndex = cidLeft; 
      return;
    }
   
    temp =  _dtlclvs[cidRight]._survivingTransferSums * (_PT[e] * freq);
    temp *= _dtlclvs[cidLeft]._uq[e];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      sampleTransferEvent(cidRight, 
          _dtlclvs[cidRight]._survivingTransferSums,
          recCell->event.pllDestSpeciesNode);
      recCell->event.destSpeciesNode = 
        recCell->event.pllDestSpeciesNode->node_index;
      recCell->event.leftGeneIndex = cidLeft; 
      recCell->event.rightGeneIndex = cidRight; 
      return;
    }
   
    // add this split contribution to the total
  }
  if (not isSpeciesLeaf) {
    // SL event
    temp = _dtlclvs[cid]._uq[f] * (_uE[g] * _PS[e]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = f;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      return;
    }
    temp = _dtlclvs[cid]._uq[g] * (_uE[f] * _PS[e]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      return;
    }
  }
}
  
template <class REAL>
REAL UndatedDTLMultiModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}



