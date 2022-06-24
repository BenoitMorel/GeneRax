#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include "MultiModel.hpp"

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDLMultiModel: public MultiModelTemplate<REAL> {
public: 
  UndatedDLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile);

  virtual ~UndatedDLMultiModel() {}

  virtual void setRates(const RatesVector &);
  virtual double computeLogLikelihood();
  virtual corax_rnode_t *sampleSpeciesNode();
  

private:
  
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  using DLCLV = std::vector<REAL>;
  std::vector<DLCLV> _dlclvs;

  REAL getLikelihoodFactor() const;
  virtual void recomputeSpeciesProbabilities();
  virtual void computeProbability(CID cid, 
    corax_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell<REAL> *recCell = nullptr);

  
};

template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &geneTreesFile):
  MultiModelTemplate<REAL>(speciesTree,
      geneSpeciesMapping,
      info,
      geneTreesFile),
  _PD(speciesTree.getNodesNumber(), 0.2),
  _PL(speciesTree.getNodesNumber(), 0.2),
  _PS(speciesTree.getNodesNumber(), 1.0),
  _uE(speciesTree.getNodesNumber(), 0.0)
{
  std::vector<REAL> zeros(speciesTree.getNodesNumber(), REAL());
  _dlclvs = std::vector<std::vector<REAL> >(
      this->_ccp.getCladesNumber(), zeros);
  for (unsigned int e = 0; e < this->_speciesTree.getNodesNumber(); ++e) {
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  this->onSpeciesTreeChange(nullptr);
}

template <class REAL>
void UndatedDLMultiModel<REAL>::setRates(const RatesVector &rates) 
{
  assert(rates.size() == 2);
  auto &dupRates = rates[0];
  auto &lossRates = rates[1];
  assert(this->_allSpeciesNodesCount == dupRates.size());
  assert(this->_allSpeciesNodesCount == lossRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PS.resize(this->_allSpeciesNodesCount);
  for (unsigned int e = 0; e < this->_allSpeciesNodesCount; ++e) {
    if (this->_info.noDup) {
      _PD[e] = 0.0;
    }
    auto sum = _PD[e] + _PL[e] + 1.0;
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] = 1.0 / sum;
  } 
  recomputeSpeciesProbabilities();
}

template <class REAL>
double UndatedDLMultiModel<REAL>::computeLogLikelihood()
{ 
  if (this->_ccp.skip()) {
    return 0.0;
  }
  this->beforeComputeLogLikelihood();
  std::vector<REAL> zeros(this->_speciesTree.getNodesNumber(), REAL());
  _dlclvs = std::vector<std::vector<REAL> >(
      this->_ccp.getCladesNumber(), zeros);
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    for (auto speciesNode: this->_allSpeciesNodes) {
      computeProbability(cid, 
          speciesNode, 
          _dlclvs[cid][speciesNode->node_index]);
    }
  }
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  REAL res = REAL();
  for (auto speciesNode: this->_allSpeciesNodes) {
    res += _dlclvs[rootCID][speciesNode->node_index];
  }
  // the root correction makes sure that UndatedDLMultiModel and
  // UndatedDL model are equivalent when there is one tree per
  // family: the UndatedDLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDL, so we multiply back here
  REAL rootCorrection(double(this->_ccp.getRootsNumber()));
  return log(res) - log(getLikelihoodFactor()) + log(rootCorrection);
}

template <class REAL>
void UndatedDLMultiModel<REAL>::recomputeSpeciesProbabilities()
{
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (this->getSpeciesLeft(speciesNode)) {
      c += _PS[e] * _uE[this->getSpeciesLeft(speciesNode)->node_index]  * 
        _uE[this->getSpeciesRight(speciesNode)->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    _uE[speciesNode->node_index] = proba;
  }
}

template <class REAL>
void UndatedDLMultiModel<REAL>::computeProbability(CID cid, 
    corax_rnode_t *speciesNode, 
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
  // terminal gene and species nodes
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
  
  for (const auto &cladeSplit: this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left; 
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    if (not isSpeciesLeaf) {
      // S event;
      temp = _dlclvs[cidLeft][f] * _dlclvs[cidRight][g] * (_PS[e] * freq); 
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        return;
      }
      temp = _dlclvs[cidRight][f] * _dlclvs[cidLeft][g] * (_PS[e] * freq); 
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        return;
      }
    }
    // D events
    temp = _dlclvs[cidLeft][e] * _dlclvs[cidRight][e] * (_PD[e] * freq);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_D;
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      return;
    }
  }
  if (not isSpeciesLeaf) {
    // SL event
    temp = _dlclvs[cid][f] * (_uE[g] * _PS[e]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = f;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      return;
    }
    
    temp = _dlclvs[cid][g] * (_uE[f] * _PS[e]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      return;
    }
  }
  // DL event
  //proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
  if (recCell) {
    std::cerr << "cerr " << proba << " " << maxProba << " " << (proba > maxProba) << std::endl;
    // we haven't sampled any event...
    assert(false);
  }
}
  
template <class REAL>
REAL UndatedDLMultiModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}


template <class REAL>
corax_rnode_t *UndatedDLMultiModel<REAL>::sampleSpeciesNode()
{
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto &uq = _dlclvs[rootCID];
  auto totalLL = std::accumulate(uq.begin(), uq.end(), REAL());
  REAL toSample = totalLL * Random::getProba();
  REAL sumLL = REAL();
  for (auto speciesNode: this->_allSpeciesNodes) {
     sumLL += uq[speciesNode->node_index];
     if (sumLL > toSample) {
        return speciesNode;
     }
  }
  assert(false);
  return nullptr;
}


 

