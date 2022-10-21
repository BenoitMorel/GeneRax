#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include "MultiModel.hpp"

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDTLMultiModel: public MultiModelTemplate<REAL> {
public: 
  UndatedDTLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile);

  virtual ~UndatedDTLMultiModel() {}

  virtual void setRates(const RatesVector &);
  virtual double computeLogLikelihood();
  virtual corax_rnode_t *sampleSpeciesNode();
  
private:
  
  
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PT; // Transfer probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  double _transferExtinctionSum;
  
  TransferConstaint _transferConstraint;

  /**
   *  All intermediate results needed to compute the reconciliation likelihood
   *  each gene node has one DTLCLV object
   *  Each DTLCLV gene  object is a function of the DTLCLVs of the direct children genes
   */
  struct DTLCLV {
    DTLCLV():
      _survivingTransferSum(REAL())
    {}

    DTLCLV(unsigned int speciesNumber):
      _uq(speciesNumber, REAL()),
      _correctionSum(speciesNumber, REAL()),
      _survivingTransferSum(REAL())
    {
    }
    // probability of a gene node rooted at a species node
    std::vector<REAL> _uq;
    std::vector<REAL> _correctionSum;

    // sum of transfer probabilities. Can be computed only once
    // for all species, to reduce computation complexity
    REAL _survivingTransferSum;
  };
  
  std::vector<DTLCLV> _dtlclvs;
  std::vector<corax_rnode_s *> _orderedSpeciations;
  std::vector<unsigned int> _orderedSpeciesRanks;



  void updateCLV(CID cid);
  REAL getLikelihoodFactor() const;
  virtual void recomputeSpeciesProbabilities();
  virtual void computeProbability(CID cid, 
    corax_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell<REAL> *recCell = nullptr);
  void sampleTransferEvent(unsigned int cid,
    corax_rnode_t *originSpeciesNode,
    Scenario::Event &event);

  
};


template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &geneTreesFile):
  MultiModelTemplate<REAL>(speciesTree,
      geneSpeciesMapping,
      info,
      geneTreesFile),
  _PD(speciesTree.getNodesNumber(), 0.2),
  _PL(speciesTree.getNodesNumber(), 0.2),
  _PT(speciesTree.getNodesNumber(), 0.2),
  _PS(speciesTree.getNodesNumber(), 1.0),
  _uE(speciesTree.getNodesNumber(), 0.0),
  _transferConstraint(info.transferConstraint)
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
  this->onSpeciesTreeChange(nullptr);
}
 
template <class REAL>
void UndatedDTLMultiModel<REAL>::setRates(const RatesVector &rates) 
{
  assert(rates.size() == 3);
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
void UndatedDTLMultiModel<REAL>::updateCLV(CID cid)
{
  auto &clv = _dtlclvs[cid];
  auto &uq = clv._uq;
  auto &correctionSum = clv._correctionSum;
  clv._survivingTransferSum = REAL();
  std::fill(uq.begin(), uq.end(), REAL());
  std::fill(correctionSum.begin(), correctionSum.end(), REAL());
  auto N = static_cast<double>(this->_allSpeciesNodes.size());
  auto sum = REAL();
  for (auto speciesNode: this->_allSpeciesNodes) {
    computeProbability(cid, 
        speciesNode, 
        uq[speciesNode->node_index]);
    sum += uq[speciesNode->node_index];
  }
  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto postOrder = this->_speciesTree.getPostOrderNodes();
    for (auto it = postOrder.rbegin();
        it != postOrder.rend(); ++it)  {
      auto speciesNode = *it;
      auto e = speciesNode->node_index;
      auto parent = speciesNode;
      while (parent) {
        auto p = parent->node_index;
        correctionSum[e] += uq[p];
        parent = parent->parent;
      }
      correctionSum[e] /= N;
    }
  }
  if (_transferConstraint == TransferConstaint::SOFTDATED) {
    std::vector<REAL> softDatedSums(N, REAL());
    REAL softDatedSum = REAL();
    for (auto leaf: this->_speciesTree.getLeaves()) {
      auto e = leaf->node_index;
      softDatedSum += uq[e];
    }
    for (auto it = this->_orderedSpeciations.rbegin(); 
        it != this->_orderedSpeciations.rend(); ++it) {
      auto node = (*it);
      auto e = node->node_index;
      softDatedSums[e] = softDatedSum;
      softDatedSum += uq[e];
    }
    for (auto node: this->_allSpeciesNodes) {
      auto e = node->node_index;
      auto p = node->parent ? node->parent->node_index : e;
      if (e != p) {
        correctionSum[e] = softDatedSums[p];
      }
      correctionSum[e] /= N;
    }
  }
  sum /= N;
  clv._survivingTransferSum = sum;
}


template <class REAL>
double UndatedDTLMultiModel<REAL>::computeLogLikelihood()
{ 
  if (this->_ccp.skip()) {
    return 0.0;
  }
  this->beforeComputeLogLikelihood();
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    updateCLV(cid);
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
    corax_rnode_t *originSpeciesNode,
    Scenario::Event &event)
{
  auto e = originSpeciesNode->node_index;
  auto &clv = _dtlclvs[cid];
  auto &survivingTransferSum = clv._survivingTransferSum;
  auto &correctionSum = clv._correctionSum;
  REAL max = REAL();
  switch (_transferConstraint) {
  case TransferConstaint::NONE:
    max = survivingTransferSum;
    break;
  case TransferConstaint::PARENTS:
    max = survivingTransferSum - correctionSum[e];
    break;
  case TransferConstaint::SOFTDATED:
    max = correctionSum[e];
    break;
  default:
    assert(false);
  }
  auto samplingProba = Random::getProba();
  max *= samplingProba;
  max *= this->_allSpeciesNodes.size();
  REAL sum = REAL();
  
  std::unordered_set<unsigned int> parents;
  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto parent = originSpeciesNode;
    while (parent) {
      parents.insert(parent->node_index);
      parent = parent->parent;
    }
  }

  for (auto speciesNode: this->_allSpeciesNodes) {
    auto h = speciesNode->node_index;
    if (_transferConstraint == TransferConstaint::NONE) {
      if (h == e) {
        continue;
      } 
    }
    // parent mode: do not continue if the receiving species
    // is a parent of the source species
    if (_transferConstraint == TransferConstaint::PARENTS) {
        if (parents.end() != parents.find(h)) {
          continue;
        }
    }
    if (_transferConstraint == TransferConstaint::SOFTDATED) {
      if (originSpeciesNode->parent) {
        auto p = originSpeciesNode->parent->node_index;
        if (_orderedSpeciesRanks[p] >= _orderedSpeciesRanks[h]) {
          continue;
        }
      }
    }
    sum += _dtlclvs[cid]._uq[h];
    if (sum > max) {
      event.pllDestSpeciesNode = speciesNode;;
      event.destSpeciesNode = h;
      return;
    }
  }
  assert(false);
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::recomputeSpeciesProbabilities()
{
  if (_transferConstraint == TransferConstaint::SOFTDATED) {
    _orderedSpeciations = this->_speciesTree.getOrderedSpeciations();  
    _orderedSpeciesRanks.resize(this->_speciesTree.getNodesNumber());
    unsigned int rank = 0;
    for (auto species: _orderedSpeciations) {
      _orderedSpeciesRanks[species->node_index] = rank++; 
    }
    for (auto leaf: this->_speciesTree.getLeaves()) {
      _orderedSpeciesRanks[leaf->node_index] = rank;
    }
  }
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
    switch (_transferConstraint) {
    case TransferConstaint::NONE:
      temp =  _dtlclvs[cidLeft]._survivingTransferSum * (_PT[e] * freq);
      break;
    case TransferConstaint::PARENTS:
      temp = (_dtlclvs[cidLeft]._survivingTransferSum - _dtlclvs[cidLeft]._correctionSum[e]) * (_PT[e] * freq);
      break;
    case TransferConstaint::SOFTDATED:
      temp = _dtlclvs[cidLeft]._correctionSum[e] * (_PT[e] * freq);
      break;
    default:
      assert(false);
    }
    temp *= _dtlclvs[cidRight]._uq[e];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      sampleTransferEvent(cidLeft, 
          speciesNode,
          recCell->event);
      recCell->event.destSpeciesNode = 
        recCell->event.pllDestSpeciesNode->node_index;
      recCell->event.leftGeneIndex = cidRight; 
      recCell->event.rightGeneIndex = cidLeft; 
      return;
    }
   
    switch (_transferConstraint) {
    case TransferConstaint::NONE:
      temp =  _dtlclvs[cidRight]._survivingTransferSum * (_PT[e] * freq);
      break;
    case TransferConstaint::PARENTS:
      temp = (_dtlclvs[cidRight]._survivingTransferSum - _dtlclvs[cidRight]._correctionSum[e]) * (_PT[e] * freq);
      break;
    case TransferConstaint::SOFTDATED:
      temp = _dtlclvs[cidRight]._correctionSum[e] * (_PT[e] * freq);
      break;
    default:
      assert(false);
    }
    temp *= _dtlclvs[cidLeft]._uq[e];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      sampleTransferEvent(cidRight, 
          speciesNode,
          recCell->event);
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

template <class REAL>
corax_rnode_t *UndatedDTLMultiModel<REAL>::sampleSpeciesNode()
{
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto &uq = _dtlclvs[rootCID]._uq;
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

