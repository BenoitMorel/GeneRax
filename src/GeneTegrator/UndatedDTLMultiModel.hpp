#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/DatedTree.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include "MultiModel.hpp"

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDTLMultiModel: public MultiModelTemplate<REAL> {
public: 
  UndatedDTLMultiModel(DatedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile);

  virtual ~UndatedDTLMultiModel() {}

  virtual void setRates(const RatesVector &);
  virtual void setAlpha(double alpha); 
  virtual double computeLogLikelihood();
  virtual corax_rnode_t *sampleSpeciesNode();
  
private:
  
  DatedTree &_datedTree;
  size_t _gammaCatNumber;
  std::vector<double> _gammaScalers;
  RatesVector _dtlRates;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PT; // Transfer probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  std::vector<double> _transferExtinctionSum;
  
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

    DTLCLV(size_t speciesNumber, size_t gammaCategories):
      _uq(speciesNumber * gammaCategories, REAL()),
      _correctionSum(speciesNumber * gammaCategories, REAL()),
      _survivingTransferSum(gammaCategories, REAL())
    {
    }
    // probability of a gene node rooted at a species node
    std::vector<REAL> _uq;
    std::vector<REAL> _correctionSum;

    // sum of transfer probabilities. Can be computed only once
    // for all species, to reduce computation complexity
    std::vector<REAL> _survivingTransferSum;
  };
  
  std::vector<DTLCLV> _dtlclvs;



  void updateCLV(CID cid);
  REAL getLikelihoodFactor() const;
  virtual void recomputeSpeciesProbabilities();
  virtual void computeProbability(CID cid, 
    corax_rnode_t *speciesNode, 
    size_t category,
    REAL &proba,
    ReconciliationCell<REAL> *recCell = nullptr);
  void sampleTransferEvent(unsigned int cid,
    corax_rnode_t *originSpeciesNode,
    size_t category,
    Scenario::Event &event);

  
};


template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(DatedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &geneTreesFile):
  MultiModelTemplate<REAL>(speciesTree.getRootedTree(),
      geneSpeciesMapping,
      info,
      geneTreesFile),
  _datedTree(speciesTree),
  _gammaCatNumber(info.gammaCategories),
  _gammaScalers(_gammaCatNumber, 1.0),
  _PD(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.2),
  _PL(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.2),
  _PT(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.2),
  _PS(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 1.0),
  _uE(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.0),
  _transferConstraint(info.transferConstraint)
{
  std::vector<REAL> zeros(this->_speciesTree.getNodesNumber());
  DTLCLV nullCLV(this->_allSpeciesNodesCount, _gammaCatNumber);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_ccp.getCladesNumber()), nullCLV);
  auto N = this->_speciesTree.getNodesNumber();
  for (size_t i = 0; i < _gammaCatNumber; ++i) {
    _gammaScalers[i] = 1.0 + 0.5 * double(i);//pow(2.0, double(i));
  }
  _gammaScalers[0] = 10;
  _dtlRates.resize(3);
  for (unsigned int i = 0; i < 3; ++i) {
    _dtlRates[i] = std::vector<double>(N, 0.2);
  }
  setAlpha(100.0);
  this->onSpeciesTreeChange(nullptr);
}
 
template <class REAL>
void UndatedDTLMultiModel<REAL>::setRates(const RatesVector &rates) 
{
  assert(rates.size() == 3);
  _dtlRates = rates;
  recomputeSpeciesProbabilities();
}
template <class REAL>
void UndatedDTLMultiModel<REAL>::setAlpha(double alpha)
{
  corax_compute_gamma_cats(alpha, _gammaScalers.size(), &_gammaScalers[0], 
      CORAX_GAMMA_RATES_MEAN);
  recomputeSpeciesProbabilities(); 
}


template <class REAL>
void UndatedDTLMultiModel<REAL>::updateCLV(CID cid)
{
  auto &clv = _dtlclvs[cid];
  auto &uq = clv._uq;
  auto &correctionSum = clv._correctionSum;
  
  std::fill(clv._survivingTransferSum.begin(), 
      clv._survivingTransferSum.end(), 
      REAL());
  std::fill(uq.begin(), uq.end(), REAL());
  std::fill(correctionSum.begin(), correctionSum.end(), REAL());
  auto N = this->_allSpeciesNodesCount;
  std::vector<REAL> sums(_gammaCatNumber, REAL());
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      auto ec = e * _gammaCatNumber + c;  
      computeProbability(cid, 
          speciesNode,
          c,
          uq[ec]);
        sums[c] += uq[ec];
    }
  }
  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto postOrder = this->_speciesTree.getPostOrderNodes();
    for (auto it = postOrder.rbegin();
        it != postOrder.rend(); ++it)  {
      auto speciesNode = *it;
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto parent = speciesNode;
        auto ec = e * _gammaCatNumber + c;
        while (parent) {
          auto p = parent->node_index;
          auto pc = p * _gammaCatNumber + c;
          correctionSum[ec] += uq[pc];
          parent = parent->parent;
        }
        correctionSum[ec] /= static_cast<double>(N);
      }
    }
  }
  if (_transferConstraint == TransferConstaint::SOFTDATED) {
    std::vector<REAL> softDatedSums(N * _gammaCatNumber, REAL());
    std::vector<REAL> softDatedSum(_gammaCatNumber, REAL());
    for (auto leaf: this->_speciesTree.getLeaves()) {
      auto e = leaf->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        softDatedSum[c] += uq[ec];
      }
    }
    for (auto it = _datedTree.getOrderedSpeciations().rbegin(); 
        it != _datedTree.getOrderedSpeciations().rend(); ++it) {
      auto node = (*it);
      auto e = node->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        softDatedSums[ec] = softDatedSum[c];
        softDatedSum[c] += uq[ec];
      }
    }
    for (auto node: this->_allSpeciesNodes) {
      auto e = node->node_index;
      auto p = node->parent ? node->parent->node_index : e;
      if (e != p) {
        for (size_t c = 0; c < _gammaCatNumber; ++c) {
          auto ec = e * _gammaCatNumber + c;
          correctionSum[ec] = softDatedSums[ec];
        }
      }
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        correctionSum[ec] /= static_cast<double>(N);
      }
    }
  }
  for (size_t c = 0; c < _gammaCatNumber; ++c) {
    sums[c] /= static_cast<double>(N);
    clv._survivingTransferSum[c] = sums[c];
  }
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
  /*
  for (auto speciesNode: this->_allSpeciesNodes) {
    res += _dtlclvs[rootCID]._uq[speciesNode->node_index];
  }
  */
  for (size_t c = 0; c < _gammaCatNumber; ++c) {
    res += _dtlclvs[rootCID]._uq[this->_speciesTree.getRoot()->node_index * _gammaCatNumber + c];
  }
  res /= double(_gammaCatNumber);
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
    size_t category,
    Scenario::Event &event)
{
  auto e = originSpeciesNode->node_index;
  auto N = this->_allSpeciesNodesCount;
  auto c = category;
  auto ec = e * _gammaCatNumber + c;
  auto &clv = _dtlclvs[cid];
  auto &survivingTransferSum = clv._survivingTransferSum;
  auto &correctionSum = clv._correctionSum;
  REAL max = REAL();
  switch (_transferConstraint) {
  case TransferConstaint::NONE:
    max = survivingTransferSum[c];
    break;
  case TransferConstaint::PARENTS:
    max = survivingTransferSum[c] - correctionSum[ec];
    break;
  case TransferConstaint::SOFTDATED:
    max = correctionSum[ec];
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
        if (_datedTree.getRank(p) >= _datedTree.getRank(h)) {
          continue;
        }
      }
    }
    auto hc = h * _gammaCatNumber + c;
    sum += _dtlclvs[cid]._uq[hc];
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
  auto &dupRates = _dtlRates[0];
  auto &lossRates = _dtlRates[1];
  auto &transferRates = _dtlRates[2];
  auto N = this->_allSpeciesNodesCount;
  assert(N == dupRates.size());
  assert(N == lossRates.size());
  assert(N == transferRates.size());
  for (unsigned int e = 0; e < this->_allSpeciesNodesCount; ++e) {
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      auto ec = e * _gammaCatNumber + c;
      _PD[ec] = dupRates[e];
      _PL[ec] = lossRates[e];
      _PT[ec] = transferRates[e];
      _PS[ec] = _gammaScalers[c];
      if (this->_info.noDup) {
        _PD[ec] = 0.0;
      }
      auto sum = _PD[ec] + _PL[ec] + _PT[ec] + _PS[ec];
      _PD[ec] /= sum;
      _PL[ec] /= sum;
      _PT[ec] /= sum;
      _PS[ec] /= sum;
    }
  }
  
  std::fill(_uE.begin(), _uE.end(), 0.0);
  
  _transferExtinctionSum = std::vector<double>(_gammaCatNumber, REAL());
  for (unsigned int it = 0; it < 4; ++it) {
    for (auto speciesNode: this->_allSpeciesNodes) {
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        double proba(_PL[ec]);
        proba += _uE[ec] * _uE[ec] * _PD[ec];
        proba += _transferExtinctionSum[c] * _PT[ec] * _uE[ec];
        if (this->getSpeciesLeft(speciesNode)) {
          auto g = this->getSpeciesLeft(speciesNode)->node_index;
          auto h = this->getSpeciesRight(speciesNode)->node_index;
          auto gc = g * _gammaCatNumber + c;
          auto hc = h * _gammaCatNumber + c;
          proba += _uE[gc]  * _uE[hc] * _PS[ec];
        }
        _uE[ec] = proba;
      }
    }
    std::fill(_transferExtinctionSum.begin(),
        _transferExtinctionSum.end(),
        0.0);
    for (auto speciesNode: this->_allSpeciesNodes) {
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        _transferExtinctionSum[c] += _uE[ec];
      }
    }
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      _transferExtinctionSum[c] /= double(N);
    }
  }
}


template <class REAL>
void UndatedDTLMultiModel<REAL>::computeProbability(CID cid, 
    corax_rnode_t *speciesNode,
    size_t category,
    REAL &proba,
    ReconciliationCell<REAL> *recCell
    )
{
  proba = REAL();
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto e = speciesNode->node_index;
  auto c = category;
  auto N = this->_allSpeciesNodesCount;
  auto ec = e * _gammaCatNumber + c;
  REAL maxProba = REAL();
  if (recCell) {
    recCell->event.geneNode = cid; 
    recCell->event.speciesNode = e;
    recCell->event.type = ReconciliationEventType::EVENT_None; 
    maxProba = recCell->maxProba;
  }
  if (this->_ccp.isLeaf(cid) && isSpeciesLeaf) {
    if (this->_geneToSpecies[cid] == e) {
      proba = REAL(_PS[ec]);
    }
    return;
  }
  REAL temp;
  unsigned int f = 0;
  unsigned int g = 0;
  unsigned int fc = 0;
  unsigned int gc = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
    fc = f * _gammaCatNumber + c;
    gc = g * _gammaCatNumber + c;
  }
 
  // for internal gene nodes
  for (const auto &cladeSplit: this->_ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left; 
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    if (not isSpeciesLeaf) {
      // S event;
      temp = _dtlclvs[cidLeft]._uq[fc] * _dtlclvs[cidRight]._uq[gc] * (_PS[ec] * freq); 
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        return;
      }
      temp = _dtlclvs[cidRight]._uq[fc] * _dtlclvs[cidLeft]._uq[gc] * (_PS[ec] * freq); 
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
    temp = _dtlclvs[cidLeft]._uq[ec] * _dtlclvs[cidRight]._uq[ec] * (_PD[ec] * freq);
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
      temp =  _dtlclvs[cidLeft]._survivingTransferSum[c] * (_PT[ec] * freq);
      break;
    case TransferConstaint::PARENTS:
      temp = (_dtlclvs[cidLeft]._survivingTransferSum[c] - _dtlclvs[cidLeft]._correctionSum[ec]) * (_PT[ec] * freq);
      break;
    case TransferConstaint::SOFTDATED:
      temp = _dtlclvs[cidLeft]._correctionSum[ec] * (_PT[ec] * freq);
      break;
    default:
      assert(false);
    }
    temp *= _dtlclvs[cidRight]._uq[ec];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      sampleTransferEvent(cidLeft, 
          speciesNode,
          c,
          recCell->event);
      recCell->event.destSpeciesNode = 
        recCell->event.pllDestSpeciesNode->node_index;
      recCell->event.leftGeneIndex = cidRight; 
      recCell->event.rightGeneIndex = cidLeft; 
      return;
    }
   
    switch (_transferConstraint) {
    case TransferConstaint::NONE:
      temp = _dtlclvs[cidRight]._survivingTransferSum[c] * (_PT[ec] * freq);
      break;
    case TransferConstaint::PARENTS:
      temp = (_dtlclvs[cidRight]._survivingTransferSum[c] 
          - _dtlclvs[cidRight]._correctionSum[ec]) * (_PT[ec] * freq);
      break;
    case TransferConstaint::SOFTDATED:
      temp = _dtlclvs[cidRight]._correctionSum[ec] * (_PT[ec] * freq);
      break;
    default:
      assert(false);
    }
    temp *= _dtlclvs[cidLeft]._uq[ec];
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_T;
      sampleTransferEvent(cidRight, 
          speciesNode,
          c,
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
    temp = _dtlclvs[cid]._uq[fc] * (_uE[gc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = f;
      recCell->event.pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      return;
    }
    temp = _dtlclvs[cid]._uq[gc] * (_uE[fc] * _PS[ec]);
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
  REAL factor(1.0);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}

template <class REAL>
corax_rnode_t *UndatedDTLMultiModel<REAL>::sampleSpeciesNode()
{
  return this->_speciesTree.getRoot();
  /*
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
  */
  assert(false);
  return nullptr;
}

