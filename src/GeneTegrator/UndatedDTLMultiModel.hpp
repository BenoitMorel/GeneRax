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
      const std::string &ccpFile);

  virtual ~UndatedDTLMultiModel() {}


  virtual void setRates(const RatesVector &);
  virtual void setAlpha(double alpha); 
  virtual double computeLogLikelihood();
  virtual corax_rnode_t *sampleSpeciesNode(unsigned int &category);
  virtual void setHighways(const std::vector<Highway> &highways) {
    _highways = highways;
    _highwaysSizes = std::vector<size_t>(this->_speciesTree.getNodesNumber(), 0);
    for (const auto &highway: _highways) {
      _highwaysSizes[highway.src->node_index]++;
    }
    resetCache();
    recomputeSpeciesProbabilities();
  }
  
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
  std::unordered_map<size_t, double> _llCache;

  std::vector<Highway> _highways;
  std::vector<size_t> _highwaysSizes;
  TransferConstaint _transferConstraint;
  OriginationStrategy _originationStrategy;
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
  double getLikelihoodFactor(unsigned int category);
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
  corax_rnode_t *getSpeciesLCA();
  size_t getHash();
  void resetCache() {_llCache.clear();}

  double getHighwayRatio() const {return 10.0;}
  double getTransferWeightNorm(unsigned int e) const {
    return _highwaysSizes[e] * getHighwayRatio() 
      + double(this->getPrunedSpeciesNodeNumber());}
};


template <class REAL>
UndatedDTLMultiModel<REAL>::UndatedDTLMultiModel(DatedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &ccpFile):
  MultiModelTemplate<REAL>(speciesTree.getRootedTree(),
      geneSpeciesMapping,
      info,
      ccpFile),
  _datedTree(speciesTree),
  _gammaCatNumber(info.gammaCategories),
  _gammaScalers(_gammaCatNumber, 1.0),
  _PD(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.2),
  _PL(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.2),
  _PT(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.1),
  _PS(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 1.0),
  _uE(this->_speciesTree.getNodesNumber() * _gammaCatNumber, 0.0),
  _transferConstraint(info.transferConstraint),
  _originationStrategy(info.originationStrategy)
{
  std::vector<REAL> zeros(this->_speciesTree.getNodesNumber());
  DTLCLV nullCLV(this->getPrunedSpeciesNodeNumber(), _gammaCatNumber);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_ccp.getCladesNumber()), nullCLV);
  auto N = this->_speciesTree.getNodesNumber();
  _highwaysSizes = std::vector<size_t>(N, 0); 
  _dtlRates.resize(3);
  for (unsigned int i = 0; i < 3; ++i) {
    _dtlRates[i] = std::vector<double>(N, 0.2);
  }
  setAlpha(1.0);
  this->onSpeciesTreeChange(nullptr);
}
 
template <class REAL>
void UndatedDTLMultiModel<REAL>::setRates(const RatesVector &rates) 
{
  resetCache();
  assert(rates.size() == 3);
  _dtlRates = rates;
  recomputeSpeciesProbabilities();
}
template <class REAL>
void UndatedDTLMultiModel<REAL>::setAlpha(double alpha)
{
  resetCache();
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
  auto N = this->getPrunedSpeciesNodeNumber();
  std::vector<REAL> sums(_gammaCatNumber, REAL());
  for (auto speciesNode: this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      auto ec = e * _gammaCatNumber + c;  
      computeProbability(cid, 
          speciesNode,
          c,
          uq[ec]);
      auto v = uq[ec];
      v /=  getTransferWeightNorm(e);
      scale(v);
      sums[c] += v;
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
          auto temp = uq[pc] / getTransferWeightNorm(p);
          scale(temp);
          correctionSum[ec] += temp;
          parent = parent->parent;
        }
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
        softDatedSum[c] += uq[ec] / getTransferWeightNorm(e);
      }
    }
    for (auto it = _datedTree.getOrderedSpeciations().rbegin(); 
        it != _datedTree.getOrderedSpeciations().rend(); ++it) {
      auto node = (*it);
      auto e = node->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        softDatedSums[ec] = softDatedSum[c];
        auto temp = uq[ec] / getTransferWeightNorm(e);
        scale(temp);
        softDatedSum[c] += temp; 
      }
    }
    for (auto node: this->getPrunedSpeciesNodes()) {
      auto e = node->node_index;
      auto p = node->parent ? node->parent->node_index : e;
      if (e != p) {
        for (size_t c = 0; c < _gammaCatNumber; ++c) {
          auto pc = p * _gammaCatNumber + c;
          auto ec = e * _gammaCatNumber + c;
          correctionSum[ec] = softDatedSums[pc];
        }
      }
    }
  }
  for (size_t c = 0; c < _gammaCatNumber; ++c) {
    clv._survivingTransferSum[c] = sums[c];
  }
}


template <class REAL>
double UndatedDTLMultiModel<REAL>::computeLogLikelihood()
{ 
  if (this->_ccp.skip()) {
    return 0.0;
  }
  auto hash = this->getHash();
  auto cacheIt = _llCache.find(hash);
  if (cacheIt != _llCache.end()) {
    return cacheIt->second;
  }

  this->beforeComputeLogLikelihood();
  for (CID cid = 0; cid < this->_ccp.getCladesNumber(); ++cid) {
    updateCLV(cid);
  }
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  std::vector<REAL> categoryLikelihoods(_gammaCatNumber, REAL());
  
  std::vector<corax_rnode_t *> rootNode;
  std::vector<corax_rnode_t *> *speciesNodes;
  switch (_originationStrategy) {
  case OriginationStrategy::ROOT:
    rootNode.push_back(this->_speciesTree.getRoot());
    speciesNodes = &rootNode;
    break;
  case OriginationStrategy::UNIFORM:
    speciesNodes = &(this->getPrunedSpeciesNodes());
    break;
  case OriginationStrategy::LCA:
    rootNode.push_back(getSpeciesLCA());
    speciesNodes = &rootNode;
    break;
  default:
    assert(false);
  }
  for (auto speciesNode: *speciesNodes) {
    auto e = speciesNode->node_index;
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      categoryLikelihoods[c] += _dtlclvs[rootCID]._uq[e * _gammaCatNumber + c];
    }
  }
  // condition on survival
  for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
    categoryLikelihoods[c] /= getLikelihoodFactor(c); 
  }
  // sum over the categories
  REAL res = std::accumulate(categoryLikelihoods.begin(), categoryLikelihoods.end(), REAL());
  // normalize by the number of categories
  res /= double(_gammaCatNumber);
  // ths root correction makes sure that UndatedDTLMultiModel and
  // UndatedDTL model are equivalent when there is one tree per
  // family: the UndatedDTLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDTL, so we multiply back here
  auto rootCorrection = double(this->_ccp.getRootsNumber()); 
  res *= rootCorrection; 
  auto ret = log(res);
  _llCache[hash] = ret; 
  return ret;
}

template <class REAL>
void UndatedDTLMultiModel<REAL>::sampleTransferEvent(unsigned int cid,
    corax_rnode_t *originSpeciesNode,
    size_t category,
    Scenario::Event &event)
{
  auto e = originSpeciesNode->node_index;
  auto N = this->getPrunedSpeciesNodeNumber();
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
  max *= N;
  REAL sum = REAL();
  
  std::unordered_set<unsigned int> parents;
  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto parent = originSpeciesNode;
    while (parent) {
      parents.insert(parent->node_index);
      parent = parent->parent;
    }
  }

  for (auto speciesNode: this->getAllSpeciesNodes()) {
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
  auto maxSpeciesId = this->_speciesTree.getNodesNumber();
  assert(maxSpeciesId == dupRates.size());
  assert(maxSpeciesId == lossRates.size());
  assert(maxSpeciesId == transferRates.size());
  for (unsigned int e = 0; e < this->getPrunedSpeciesNodeNumber(); ++e) {
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
  auto transferSum = std::vector<double>(_gammaCatNumber, REAL());
  for (unsigned int it = 0; it < 4; ++it) {
    for (auto speciesNode: this->getAllSpeciesNodes()) {
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        double proba(_PL[ec]);
        proba += _uE[ec] * _uE[ec] * _PD[ec];
        proba += transferSum[c] * _PT[ec] * _uE[ec];
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
    std::fill(transferSum.begin(),
        transferSum.end(),
        0.0);
    for (auto speciesNode: this->getAllSpeciesNodes()) {
      auto e = speciesNode->node_index;
      for (size_t c = 0; c < _gammaCatNumber; ++c) {
        auto ec = e * _gammaCatNumber + c;
        transferSum[c] += _uE[ec];
      }
    }
    for (size_t c = 0; c < _gammaCatNumber; ++c) {
      transferSum[c] /= double(this->getPrunedSpeciesNodeNumber());
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
    
    // highway transfers
    for (auto highway: _highways) {
      if (e == highway.src->node_index) {
        // e is the src species, d is the receiving species
        auto d = highway.dest->node_index;
        auto dc = d * _gammaCatNumber + c;  
        temp = (_dtlclvs[cidLeft]._uq[ec] * _dtlclvs[cidRight]._uq[dc]) * (_PT[ec] * freq * getHighwayRatio() / getTransferWeightNorm(e));
        scale(temp);
        proba += temp;
        if (recCell && proba > maxProba) {
          recCell->event.type = ReconciliationEventType::EVENT_T;
          recCell->event.destSpeciesNode = d;
          recCell->event.leftGeneIndex = cidLeft;
          recCell->event.rightGeneIndex = cidRight; 
        }
        temp = (_dtlclvs[cidRight]._uq[ec] * _dtlclvs[cidLeft]._uq[dc]) * (_PT[ec] * freq * getHighwayRatio() / getTransferWeightNorm(e));
        scale(temp);
        proba += temp;
        if (recCell && proba > maxProba) {
          recCell->event.type = ReconciliationEventType::EVENT_T;
          recCell->event.destSpeciesNode = d;
          recCell->event.leftGeneIndex = cidRight;
          recCell->event.rightGeneIndex = cidLeft; 
        }
      }
    }
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
      recCell->event.pllLostSpeciesNode = this->getSpeciesRight(speciesNode);
      return;
    }
    temp = _dtlclvs[cid]._uq[gc] * (_uE[fc] * _PS[ec]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = g;
      recCell->event.pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      recCell->event.pllLostSpeciesNode = this->getSpeciesLeft(speciesNode);
      return;
    }
  }
}
  
template <class REAL>
double UndatedDTLMultiModel<REAL>::getLikelihoodFactor(unsigned int category) 
{
  double factor(0.0);
  for (auto speciesNode: this->getPrunedSpeciesNodes()) {
    auto e = speciesNode->node_index;
    factor += 1.0 - _uE[e * _gammaCatNumber + category];
  }
  return factor;
}

template <class REAL>
corax_rnode_t *UndatedDTLMultiModel<REAL>::sampleSpeciesNode(unsigned int &category)
{
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto &uq = _dtlclvs[rootCID]._uq;
  std::vector<corax_rnode_t *> rootNodeVector;
  std::vector<corax_rnode_t *> *speciesNodeCandidates = nullptr;
  switch (_originationStrategy) {
  case OriginationStrategy::UNIFORM:
    speciesNodeCandidates = &(this->getPrunedSpeciesNodes());
    break;
  case OriginationStrategy::ROOT:
    rootNodeVector.push_back(this->_speciesTree.getRoot());
    speciesNodeCandidates = &rootNodeVector;
    break;
  case OriginationStrategy::LCA:
    rootNodeVector.push_back(getSpeciesLCA());
    speciesNodeCandidates = &rootNodeVector;
    break;
  default:
    assert(false);
  }
  assert(speciesNodeCandidates);
  REAL totalLL = REAL();
  for (auto node: *speciesNodeCandidates) {
    auto e = node->node_index;
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      totalLL += uq[e * _gammaCatNumber + c];
    }
  }
  auto toSample = totalLL * Random::getProba();
  auto sumLL = REAL();
  for (auto node: *speciesNodeCandidates) {
    auto e = node->node_index;
    for (unsigned int c = 0; c < _gammaCatNumber; ++c) {
      sumLL += uq[e * _gammaCatNumber + c];
      if (sumLL >= toSample) {
        category = c;
        return node;
      }
    }
  } 
  assert(false);
  return nullptr;
}
template <class REAL>  
corax_rnode_t *UndatedDTLMultiModel<REAL>::getSpeciesLCA()
{
  auto tree = &this->_speciesTree;
  corax_rnode_t *lca = nullptr;
  for (auto it: this->_speciesNameToId) {
    auto spid = it.second;
    auto node = tree->getNode(spid);
    lca = tree->getLCA(lca, node);
  }
  return lca;
}
  
template <class REAL> 
size_t UndatedDTLMultiModel<REAL>::getHash()
{
  auto hash = this->getSpeciesTreeHash();
  switch (_transferConstraint) {
  case TransferConstaint::NONE:
  case TransferConstaint::PARENTS:
    return hash;
  case TransferConstaint::SOFTDATED:
    return this->_datedTree.getOrderingHash(hash);
  default:
    assert(false);
  }
  return hash;
}

