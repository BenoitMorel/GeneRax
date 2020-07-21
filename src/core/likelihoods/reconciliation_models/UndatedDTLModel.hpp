#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>




/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In addition, we forbid transfers to parent species
*/
template <class REAL>
class UndatedDTLModel: public AbstractReconciliationModel<REAL> {
public:
  UndatedDTLModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMappingp, 
      bool rootedGeneTree, 
      double minGeneBranchLength,
      bool pruneSpeciesTree):
    
    AbstractReconciliationModel<REAL>(speciesTree, geneSpeciesMappingp, rootedGeneTree, minGeneBranchLength, pruneSpeciesTree)
  {
  } 
  UndatedDTLModel(const UndatedDTLModel &) = delete;
  UndatedDTLModel & operator = (const UndatedDTLModel &) = delete;
  UndatedDTLModel(UndatedDTLModel &&) = delete;
  UndatedDTLModel & operator = (UndatedDTLModel &&) = delete;
  virtual ~UndatedDTLModel();
  
  // overloaded from parent
  virtual void setRates(const RatesVector &rates);
  
  virtual void rollbackToLastState();
protected:
  // overloaded from parent
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  // overloaded from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overload from parent
  virtual void recomputeSpeciesProbabilities();
  // overloaded from parent
  virtual REAL getGeneRootLikelihood(pll_unode_t *root) const;
  // overload from parent
  virtual void computeGeneRootLikelihood(pll_unode_t *virtualRoot);
  virtual REAL getGeneRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {
    return _dtlclvs[root->node_index + this->_maxGeneId + 1]._uq[speciesRoot->node_index];
  }
  virtual REAL getLikelihoodFactor() const;
  virtual void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr,
      bool stochastic = false);
private:
  // model
  std::vector<double> _PD; // Duplication probability, per branch
  std::vector<double> _PL; // Loss probability, per branch
  std::vector<double> _PT; // Transfer probability, per branch
  std::vector<double> _PS; // Speciation probability, per branch
  // SPECIES
  std::vector<REAL> _uE; // Probability for a gene to become extinct on each brance
  REAL _transferExtinctionSum;

  
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
    {}
    // probability of a gene node rooted at a species node
    std::vector<REAL> _uq;
    // sum of transfer probabilities. Can be computed only once
    // for all species, to reduce computation complexity
    REAL _survivingTransferSums;
  };

  // Current DTLCLV values
  std::vector<DTLCLV> _dtlclvs;
private:
  void updateTransferSums(REAL &transferExtinctionSum,
      const std::vector<REAL> &probabilities);
  void getBestTransfer(pll_unode_t *parentGeneNode, 
    pll_rnode_t *originSpeciesNode,
    bool isVirtualRoot,
    pll_unode_t *&transferedGene,
    pll_unode_t *&stayingGene,
    pll_rnode_t *&recievingSpecies,
    REAL &proba, 
    bool stochastic = false);
  void getBestTransferLoss(Scenario &scenario,
      pll_unode_t *parentGeneNode, 
    pll_rnode_t *originSpeciesNode,
    pll_rnode_t *&recievingSpecies,
    REAL &proba,
    bool stochastic = false);
  unsigned int getIterationsNumber() const {return 4;}    
  REAL getCorrectedTransferExtinctionSum(unsigned int speciesId) const {
    return _transferExtinctionSum * _PT[speciesId];
  }

  REAL getCorrectedTransferSum(unsigned int geneId, unsigned int speciesId) const
  {
    return _dtlclvs[geneId]._survivingTransferSums * _PT[speciesId];
  }
  std::vector<pll_rnode_s *> &getSpeciesNodesToUpdate() {
    return this->_speciesNodesToUpdate;
  }
  std::vector<pll_rnode_s *> &getSpeciesNodesToUpdateSafe() {
    return this->_allSpeciesNodes;
  }
};


template <class REAL>
void UndatedDTLModel<REAL>::setInitialGeneTree(PLLUnrootedTree &tree)
{
  AbstractReconciliationModel<REAL>::setInitialGeneTree(tree);
  DTLCLV nullCLV(this->_allSpeciesNodesCount);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_maxGeneId + 1), nullCLV);
}

template <class REAL>
void UndatedDTLModel<REAL>::updateTransferSums(REAL &transferSum,
    const std::vector<REAL> &probabilities)
{
  transferSum = REAL();
  for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
    transferSum += probabilities[speciesNode->node_index];
  }
  transferSum /= this->_allSpeciesNodes.size();
}


template <class REAL>
void UndatedDTLModel<REAL>::setRates(const RatesVector &rates)
{
  this->_geneRoot = 0;
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
    auto sum = _PD[e] + _PL[e] + _PT[e] + 1.0;
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] = 1.0 / sum;
  } 
  recomputeSpeciesProbabilities();
  this->invalidateAllCLVs();
  this->invalidateAllSpeciesCLVs();
}

template <class REAL>
UndatedDTLModel<REAL>::~UndatedDTLModel() { }

template <class REAL>
void UndatedDTLModel<REAL>::recomputeSpeciesProbabilities()
{
  _uE.resize(this->_allSpeciesNodesCount);
  for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
    _uE[speciesNode->node_index] = REAL(0.0);
  }
  _transferExtinctionSum = REAL();
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
      auto e = speciesNode->node_index;
      if (it + 1 == getIterationsNumber() && !speciesNode->left) {
        _uE[e] = _uE[e] * (1.0 - this->_fm[e]) + REAL(this->_fm[e]);
        continue;
      }
      REAL proba(_PL[e]);
      REAL temp = _uE[e] * _uE[e] * _PD[e];
      scale(temp);
      proba += temp;
      temp = getCorrectedTransferExtinctionSum(e) * _uE[e];
      scale(temp);
      proba += temp;
      if (this->getSpeciesLeft(speciesNode)) {
        temp = _uE[this->getSpeciesLeft(speciesNode)->node_index]  * _uE[this->getSpeciesRight(speciesNode)->node_index] * _PS[e];
        scale(temp);
        proba += temp;
      }
      //PRINT_ERROR_PROBA(proba)
      _uE[speciesNode->node_index] = proba;
    }
    updateTransferSums(_transferExtinctionSum, _uE);
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::updateCLV(pll_unode_t *geneNode)
{
  auto gid = geneNode->node_index; 
  auto lca = this->_geneToSpeciesLCA[gid];
  auto lcaRight = lca;
  auto lcaLeft = lca;
  if (geneNode->next) {
    auto geneLeft = this->getLeft(geneNode, false);
    auto geneRight = this->getRight(geneNode, false);
    lcaRight = this->_geneToSpeciesLCA[geneLeft->node_index];
    lcaLeft = this->_geneToSpeciesLCA[geneRight->node_index];
  }
  auto &ancestorsLeft = this->_speciesTree.getAncestorssCache(lcaLeft);
  auto &ancestorsRight = this->_speciesTree.getAncestorssCache(lcaRight);
  
  _dtlclvs[gid]._survivingTransferSums = REAL();
  for (auto speciesNode: getSpeciesNodesToUpdate()) { 
    auto e = speciesNode->node_index;
    if (!(ancestorsLeft[e] || ancestorsRight[e])) {
      _dtlclvs[gid]._uq[speciesNode->node_index] = REAL();
    }
  }
  
  for (auto speciesNode: getSpeciesNodesToUpdateSafe()) { 
    auto e = speciesNode->node_index;
    //if (parentsCache[e]) {
    if (ancestorsLeft[e] || ancestorsRight[e]) {
      computeProbability(geneNode, 
        speciesNode, 
        _dtlclvs[gid]._uq[e]);
    }
  }
  updateTransferSums(_dtlclvs[gid]._survivingTransferSums, 
    _dtlclvs[gid]._uq);
}


template <class REAL>
void UndatedDTLModel<REAL>::computeGeneRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  _dtlclvs[u]._survivingTransferSums = REAL();
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    auto e = speciesNode->node_index;
    _dtlclvs[u]._uq[e] = REAL();
  }
  auto geneLeft = this->getLeft(virtualRoot, true);
  auto geneRight = this->getRight(virtualRoot, true);
  auto lcaRight = this->_geneToSpeciesLCA[geneLeft->node_index];
  auto lcaLeft = this->_geneToSpeciesLCA[geneRight->node_index];
  auto &ancestorsLeft = this->_speciesTree.getAncestorssCache(lcaLeft);
  auto &ancestorsRight = this->_speciesTree.getAncestorssCache(lcaRight);
  
  
  for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
    unsigned int e = speciesNode->node_index;
    if (ancestorsLeft[e] || ancestorsRight[e]) {
      computeProbability(virtualRoot, speciesNode, _dtlclvs[u]._uq[e], true);
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot,
      Scenario *scenario,
      Scenario::Event *event,
      bool stochastic)
{
  
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  
  if (event) {
    event->geneNode = gid; 
    event->speciesNode = e;
    event->type = ReconciliationEventType::EVENT_None; 
  }

  if (isSpeciesLeaf and isGeneLeaf and e == this->_geneToSpecies[gid]) {
    proba = REAL(_PS[e]);
    return;
  }
  typedef std::array<REAL, 7>  ValuesArray;
  ValuesArray values;
  if (event) {
    values[0] = values[1] = values[2] = values[3] = REAL();
    values[4] = values[5] = values[6] = REAL();
  }
  proba = REAL();
  
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }

  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
  }

  if (not isGeneLeaf) {
    // S event
    auto u_left = leftGeneNode->node_index;
    auto u_right = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      //  speciation event
      values[0] = _dtlclvs[u_left]._uq[f];
      values[1] = _dtlclvs[u_left]._uq[g];
      values[0] *= _dtlclvs[u_right]._uq[g];
      values[1] *= _dtlclvs[u_right]._uq[f];
      values[0] *= _PS[e]; 
      values[1] *= _PS[e]; 
      scale(values[0]);
      scale(values[1]);
      proba += values[0];
      proba += values[1];
    }
    // D event
    values[2] = _dtlclvs[u_left]._uq[e];
    values[2] *= _dtlclvs[u_right]._uq[e];
    values[2] *= _PD[e];
    scale(values[2]);
    proba += values[2];
    
    // T event
    values[5] = getCorrectedTransferSum(u_left, e);
    values[5] *= _dtlclvs[u_right]._uq[e];
    scale(values[5]);
    values[6] = getCorrectedTransferSum(u_right, e);
    values[6] *= _dtlclvs[u_left]._uq[e];
    scale(values[6]);
    proba += values[5];
    proba += values[6];
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _dtlclvs[gid]._uq[f];
    values[3] *= (_uE[g] * _PS[e]);
    scale(values[3]);
    values[4] = _dtlclvs[gid]._uq[g];
    values[4]*= _uE[f] * _PS[e];
    scale(values[4]);
    proba += values[3];
    proba += values[4];
  }
  if (event) {
    assert(scenario);
    pll_unode_t *transferedGene = 0;
    pll_unode_t *stayingGene = 0;
    pll_rnode_t *recievingSpecies = 0;
    values[5] = values[6] = REAL(); // invalidate these ones
    if (!isGeneLeaf) {
      getBestTransfer(geneNode, speciesNode, isVirtualRoot, 
          transferedGene, stayingGene, recievingSpecies, values[5], stochastic);
    }
    int maxValueIndex = 0;
    if (!stochastic) {
      maxValueIndex =static_cast<unsigned int>(std::distance(values.begin(),
          std::max_element(values.begin(), values.end())
          ));
    } else {
      maxValueIndex = sampleIndex<ValuesArray, REAL>(values);
    }
    if (-1 == maxValueIndex || values[maxValueIndex] == REAL()) {
      event->type = ReconciliationEventType::EVENT_Invalid;
      return;
    }
    switch(maxValueIndex) {
    case 0:
      event->type = ReconciliationEventType::EVENT_S;
      event->cross = false;
      break;
    case 1:
      event->type = ReconciliationEventType::EVENT_S;
      event->cross = true;
      break;
    case 2:
      event->type = ReconciliationEventType::EVENT_D;
      break;
    case 3:
      event->type = ReconciliationEventType::EVENT_SL;
      event->destSpeciesNode = f;
      event->pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      break;
    case 4:
      event->type = ReconciliationEventType::EVENT_SL;
      event->destSpeciesNode = g;
      event->pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      break;
    case 5:
      event->type = ReconciliationEventType::EVENT_T;
      event->transferedGeneNode = transferedGene->node_index;
      event->destSpeciesNode = recievingSpecies->node_index;
      event->pllTransferedGeneNode = transferedGene;
      event->pllDestSpeciesNode = recievingSpecies;
      break;
    case 6:
      assert(false);
      break;
    default:
      assert(false);
    };
  }
}



template <class REAL>
REAL UndatedDTLModel<REAL>::getGeneRootLikelihood(pll_unode_t *root) const
{
  REAL sum = REAL();
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    sum += _dtlclvs[u]._uq[e];
  }
  return sum;
}

template <class REAL>
REAL UndatedDTLModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - _uE[e]);
  }
  return factor;
      
}


template <class REAL>
void UndatedDTLModel<REAL>::rollbackToLastState()
{
}

template <class REAL>
void UndatedDTLModel<REAL>::getBestTransfer(pll_unode_t *parentGeneNode, 
  pll_rnode_t *originSpeciesNode,
  bool isVirtualRoot,
  pll_unode_t *&transferedGene,
  pll_unode_t *&stayingGene,
  pll_rnode_t *&recievingSpecies,
  REAL &proba, 
  bool stochastic)
{
  unsigned int speciesNumber = this->_allSpeciesNodes.size();
  proba = REAL();
  auto e = originSpeciesNode->node_index;
  auto u_left = this->getLeft(parentGeneNode, isVirtualRoot);
  auto u_right = this->getRight(parentGeneNode, isVirtualRoot);
  std::vector<REAL> transferProbas(speciesNumber * 2, REAL());
  double factor = _PT[e] / static_cast<double>(speciesNumber);
  for (auto species: this->_allSpeciesNodes) {
    auto h = species->node_index;
    if (h == e) {
      continue;
    }
    transferProbas[h] = (_dtlclvs[u_left->node_index]._uq[h] 
        * _dtlclvs[u_right->node_index]._uq[e]) * factor;
    transferProbas[h + speciesNumber] = (_dtlclvs[u_right->node_index]._uq[h] 
        * _dtlclvs[u_left->node_index]._uq[e]) * factor;
  }
  if (stochastic) {
    // stochastic sample: proba will be set to the sum of probabilities
    proba = REAL();
    for (auto &value: transferProbas) {
      proba += value;
    }
    auto bestIndex = sampleIndex<std::vector<REAL>, REAL>(transferProbas);
    if (bestIndex == -1) {
      proba = REAL();
      return;
    }
    bool left = static_cast<unsigned int>(bestIndex) < speciesNumber;
    transferedGene = left ? u_left : u_right;
    stayingGene = !left ? u_left : u_right;
    // I am not sure I can find the species from
    // its index in the species array, so I keep it safe here
    for (auto species: this->_allSpeciesNodes) {
      if (species->node_index == bestIndex % speciesNumber) {
        recievingSpecies = species;
      }
    }
  } else {
    // find the max
    for (auto species: this->_allSpeciesNodes) {
      auto h = species->node_index;
      if (proba < transferProbas[h]) {
        proba = transferProbas[h];
        transferedGene = u_left;
        stayingGene = u_right;
        recievingSpecies = species;
      }
      if (proba < transferProbas[h + speciesNumber]) {
        proba = transferProbas[h + speciesNumber];
        transferedGene = u_right;
        stayingGene = u_left;
        recievingSpecies = species;
      }
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::getBestTransferLoss(Scenario &scenario,
   pll_unode_t *parentGeneNode, 
  pll_rnode_t *originSpeciesNode,
  pll_rnode_t *&recievingSpecies,
  REAL &proba, 
  bool stochastic)
{
  proba = REAL();
  auto e = originSpeciesNode->node_index;
  auto u = parentGeneNode->node_index;
  
  unsigned int speciesNumber = this->_allSpeciesNodes.size();
  std::vector<REAL> transferProbas(speciesNumber, REAL());
  REAL factor = _uE[e] * (_PT[e] / static_cast<double>(this->_allSpeciesNodes.size()));
  for (auto species: this->_allSpeciesNodes) {
    auto h = species->node_index;
    if (h == e) {
      continue;
    }
    transferProbas[h] = _dtlclvs[u]._uq[h] * factor;
  }
  if (!stochastic) {
    for (auto species: this->_allSpeciesNodes) {
      auto h = species->node_index;
      if (proba < transferProbas[h]) {
        if (!scenario.isBlacklisted(u, h)) {
          scenario.blackList(u, h);
          proba = transferProbas[h];
          recievingSpecies = species;
        }
      }
    } 
  } else {
    proba = REAL();
    for (auto p: transferProbas) {
      proba += p;
    }
    unsigned int h = 0;
    do {

      int bestIndex = sampleIndex<std::vector<REAL>, REAL>(transferProbas);
      if (bestIndex == -1) {
        proba = REAL();
        return;
      }
      for (auto species: this->_allSpeciesNodes) {
        h = species->node_index;
        if (static_cast<unsigned int>(bestIndex) == h) {
          recievingSpecies = species;
          transferProbas[h] = REAL(); // in case it's blacklisted, avoid infinite loop
          break;
        }
      }
    } while (scenario.isBlacklisted(u, h));
    scenario.blackList(u, h);
  }
}
