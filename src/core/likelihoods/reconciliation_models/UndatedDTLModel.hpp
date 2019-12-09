#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>

#define IS_PROBA(x) true //(REAL(0.0) <= (x) && (x) <= REAL(1.0))
  
#define PRINT_ERROR_PROBA(x) // if (!IS_PROBA(x)) {std::cerr << "error " << x << std::endl;} assert(IS_PROBA(x));  



/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In addition, we forbid transfers to parent species
*/
template <class REAL>
class UndatedDTLModel: public AbstractReconciliationModel<REAL> {
public:
  UndatedDTLModel(PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMappingp, bool rootedGeneTree):
    AbstractReconciliationModel<REAL>(speciesTree, geneSpeciesMappingp, rootedGeneTree)
  {
  } 
  UndatedDTLModel(const UndatedDTLModel &) = delete;
  UndatedDTLModel & operator = (const UndatedDTLModel &) = delete;
  UndatedDTLModel(UndatedDTLModel &&) = delete;
  UndatedDTLModel & operator = (UndatedDTLModel &&) = delete;
  virtual ~UndatedDTLModel();
  
  // overloaded from parent
  virtual void setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates);
  
  virtual void rollbackToLastState();
protected:
  // overloaded from parent
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // overloaded from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overload from parent
  virtual void recomputeSpeciesProbabilities();
  // overloaded from parent
  virtual REAL getRootLikelihood(pll_unode_t *root) const;
  // overload from parent
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot);
  virtual REAL getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {
    return _dtlclvs[root->node_index + this->_maxGeneId + 1]._uq[speciesRoot->node_index];
  }
  virtual REAL getLikelihoodFactor() const;
  virtual void beforeComputeLogLikelihood(); 
  virtual void afterComputeLogLikelihood(); 
  virtual void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr);
private:
  // model
  std::vector<double> _PD; // Duplication probability, per branch
  std::vector<double> _PL; // Loss probability, per branch
  std::vector<double> _PT; // Transfer probability, per branch
  std::vector<double> _PS; // Speciation probability, per branch
  // SPECIES
  std::vector<REAL> _uE; // Probability for a gene to become extinct on each brance
  REAL _transferExtinctionSum;
  REAL _transferExtinctionSumBackup;

  
  /**
   *  All intermediate results needed to compute the reconciliation likelihood
   *  each gene node has one DTLCLV object
   *  Each DTLCLV gene  object is a function of the DTLCLVs of the direct children genes
   */
  struct DTLCLV {
    DTLCLV():
      _survivingTransferSums(REAL()),
      _survivingTransferSumsInvariant(REAL()),
      _survivingTransferSumsOneMore(REAL())
    {}

    DTLCLV(unsigned int speciesNumber):
      _uq(speciesNumber, REAL()),
      _survivingTransferSums(REAL()),
      _survivingTransferSumsInvariant(REAL()),
      _survivingTransferSumsOneMore(REAL())
    {}
    // probability of a gene node rooted at a species node
    std::vector<REAL> _uq;
    // sum of transfer probabilities. Can be computed only once
    // for all species, to reduce computation complexity
    REAL _survivingTransferSums;
    // subsum of transfer probbabilities that did not change
    // in case of partial likelihood recomputation
    REAL _survivingTransferSumsInvariant;
    // when computing the likelihood in slow mode, we update this value,
    // because we need it to compute _survivingTransferSumsInvariant
    // consistently in fast mode
    REAL _survivingTransferSumsOneMore;
  };

  // Current DTLCLV values
  std::vector<DTLCLV> _dtlclvs;
  // Previous DTLCVL values, to rollback to a consistent state
  // after a fast likelihood computation
  std::vector<DTLCLV> _dtlclvsBackup;
private:
  void updateTransferSums(REAL &transferExtinctionSum,
    const REAL &transferSumBackup,
    const std::vector<REAL> &probabilities);
  void resetTransferSums(const REAL &transferSum,
    REAL &transferSumBackup,
    const std::vector<REAL> &probabilities);
  void getBestTransfer(pll_unode_t *parentGeneNode, 
    pll_rnode_t *originSpeciesNode,
    bool isVirtualRoot,
    pll_unode_t *&transferedGene,
    pll_unode_t *&stayingGene,
    pll_rnode_t *&recievingSpecies,
    REAL &proba);
  void getBestTransferLoss(Scenario &scenario,
      pll_unode_t *parentGeneNode, 
    pll_rnode_t *originSpeciesNode,
    pll_rnode_t *&recievingSpecies,
    REAL &proba);
  unsigned int getIterationsNumber() const { return this->_fastMode ? 1 : 5;}    
  REAL getCorrectedTransferExtinctionSum(unsigned int speciesId) const {
    return _transferExtinctionSum * _PT[speciesId];
  }

  REAL getCorrectedTransferSum(unsigned int geneId, unsigned int speciesId) const
  {
    return _dtlclvs[geneId]._survivingTransferSums * _PT[speciesId];
  }
  std::vector<pll_rnode_s *> &getSpeciesNodesToUpdate() {
    return (this->_fastMode ? this->_speciesNodesToUpdate : this->_allSpeciesNodes);
  }
};


const unsigned int CACHE_SIZE = 100000;


template <class REAL>
void UndatedDTLModel<REAL>::setInitialGeneTree(pll_utree_t *tree)
{
  AbstractReconciliationModel<REAL>::setInitialGeneTree(tree);
  assert(this->_allSpeciesNodesCount);
  assert(this->_maxGeneId);
  std::vector<REAL> zeros(this->_allSpeciesNodesCount);
  DTLCLV nullCLV(this->_allSpeciesNodesCount);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_maxGeneId + 1), nullCLV);
  _dtlclvsBackup = std::vector<DTLCLV>(2 * (this->_maxGeneId + 1), nullCLV);
}

  template <class REAL>
void UndatedDTLModel<REAL>::resetTransferSums(const REAL &transferSum,
    REAL &transferSumInvariant,
    const std::vector<REAL> &probabilities)
{
  if (this->_fastMode) {
    REAL diff = REAL();
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      diff += probabilities[speciesNode->node_index];
    }
    diff /= this->_allSpeciesNodes.size();
    transferSumInvariant = transferSum - diff;
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::updateTransferSums(REAL &transferSum,
    const REAL &transferSumInvariant,
    const std::vector<REAL> &probabilities)
{
  transferSum = REAL();
  for (auto speciesNode:  getSpeciesNodesToUpdate()) {
    auto e = speciesNode->node_index;
    transferSum += probabilities[e];
  }
  transferSum /= this->_allSpeciesNodes.size();
  if (this->_fastMode) {
    transferSum += transferSumInvariant;
  }
}


template <class REAL>
void UndatedDTLModel<REAL>::setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates)
{
  this->_geneRoot = 0;
  assert(this->_allSpeciesNodesCount == dupRates.size());
  assert(this->_allSpeciesNodesCount == lossRates.size());
  assert(this->_allSpeciesNodesCount == transferRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PT = transferRates;
  _PS.resize(this->_allSpeciesNodesCount);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
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
  REAL unused = REAL(1.0);
  resetTransferSums(_transferExtinctionSum, unused, _uE);
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      auto e = speciesNode->node_index;
      REAL proba(_PL[e]);
      REAL temp = _uE[e] * _uE[e] * _PD[e];
      scale(temp);
      proba += temp;
      temp = getCorrectedTransferExtinctionSum(e) * _uE[e];
      scale(temp);
      proba += temp;
      if (speciesNode->left) {
        temp = _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index] * _PS[e];
        scale(temp);
        proba += temp;
      }
      //PRINT_ERROR_PROBA(proba)
      _uE[speciesNode->node_index] = proba;
    }
    updateTransferSums(_transferExtinctionSum, unused, _uE);
  }
}


template <class REAL>
void UndatedDTLModel<REAL>::updateCLV(pll_unode_t *geneNode)
{
  auto gid = geneNode->node_index;
  resetTransferSums(this->_fastMode ? _dtlclvs[gid]._survivingTransferSumsOneMore : _dtlclvs[gid]._survivingTransferSums, _dtlclvs[gid]._survivingTransferSumsInvariant, _dtlclvs[gid]._uq);
  
  if (!this->_fastMode) {
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      _dtlclvs[gid]._uq[speciesNode->node_index] = REAL();
    }
  }
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    updateTransferSums(_dtlclvs[gid]._survivingTransferSums, _dtlclvs[gid]._survivingTransferSumsInvariant, _dtlclvs[gid]._uq);
    for (auto speciesNode: getSpeciesNodesToUpdate()) { 
      computeProbability(geneNode, 
          speciesNode, 
          _dtlclvs[gid]._uq[speciesNode->node_index]);
    }
  }
  if (this->_likelihoodMode == PartialLikelihoodMode::PartialSpecies && !this->_fastMode) {
    updateTransferSums(_dtlclvs[gid]._survivingTransferSumsOneMore, _dtlclvs[gid]._survivingTransferSumsInvariant, _dtlclvs[gid]._uq);
  }
}



template <class REAL>
void UndatedDTLModel<REAL>::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot,
      Scenario *scenario,
      Scenario::Event *event)
{
  
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  

  if (event) {
    event->geneNode = gid; 
    event->speciesNode = e;
    event->type = ReconciliationEventType::EVENT_None; 
  }

  if (isSpeciesLeaf and isGeneLeaf and e == this->_geneToSpecies[gid]) {
    proba = REAL(_PS[e]);
    return;
  }
  
  std::array<REAL, 8> values;
  values[0] = values[1] = values[2] = values[3] = REAL();
  values[4] = values[5] = values[6] = values[7] = REAL();
  
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
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
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
  // TL event
  values[7] = getCorrectedTransferSum(gid, e);
  values[7] *= _uE[e];
  scale(values[7]);
  proba += values[7];
  
  if (event) {
    assert(scenario);
    pll_unode_t *transferedGene = 0;
    pll_unode_t *stayingGene = 0;
    pll_rnode_t *recievingSpecies = 0;
    pll_rnode_t *tlRecievingSpecies = 0;
    values[5] = values[6] = values[7] = REAL(); // invalidate these ones
    if (!isGeneLeaf) {
      getBestTransfer(geneNode, speciesNode, isVirtualRoot, 
          transferedGene, stayingGene, recievingSpecies, values[5]);
    }
    getBestTransferLoss(*scenario, geneNode, speciesNode, tlRecievingSpecies, values[7]);
    unsigned int maxValueIndex = static_cast<unsigned int>(std::distance(values.begin(),
          std::max_element(values.begin(), values.end())
          ));
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
      event->pllDestSpeciesNode = speciesNode->left;
      break;
    case 4:
      event->type = ReconciliationEventType::EVENT_SL;
      event->destSpeciesNode = g;
      event->pllDestSpeciesNode = speciesNode->right;
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
    case 7:
      event->type = ReconciliationEventType::EVENT_TL;
      event->transferedGeneNode = gid;
      event->destSpeciesNode = tlRecievingSpecies->node_index;
      event->pllTransferedGeneNode = geneNode;
      event->pllDestSpeciesNode = tlRecievingSpecies;
      break;
    default:
      assert(false);
    };
  }
}


template <class REAL>
void UndatedDTLModel<REAL>::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  resetTransferSums(this->_fastMode ? _dtlclvs[u]._survivingTransferSumsOneMore : _dtlclvs[u]._survivingTransferSums, _dtlclvs[u]._survivingTransferSumsInvariant, _dtlclvs[u]._uq);
  if (!this->_fastMode) {
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      auto e = speciesNode->node_index;
      _dtlclvs[u]._uq[e] = REAL();
    }
  }
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    updateTransferSums(_dtlclvs[u]._survivingTransferSums, _dtlclvs[u]._survivingTransferSumsInvariant, _dtlclvs[u]._uq);
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      unsigned int e = speciesNode->node_index;
      computeProbability(virtualRoot, speciesNode, _dtlclvs[u]._uq[e], true);
    }
  }
  if (!this->_fastMode) {
    updateTransferSums(_dtlclvs[u]._survivingTransferSumsOneMore, _dtlclvs[u]._survivingTransferSumsInvariant, _dtlclvs[u]._uq);
  }
}


template <class REAL>
REAL UndatedDTLModel<REAL>::getRootLikelihood(pll_unode_t *root) const
{
  REAL sum = REAL();
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    sum += _dtlclvs[u]._uq[e];
  }
  PRINT_ERROR_PROBA(sum);
  assert(IS_PROBA(sum));
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
void UndatedDTLModel<REAL>::beforeComputeLogLikelihood()
{
  AbstractReconciliationModel<REAL>::beforeComputeLogLikelihood();
  if (this->_likelihoodMode == PartialLikelihoodMode::PartialSpecies) {
    if (this->_fastMode) {
      _transferExtinctionSumBackup = _transferExtinctionSum;
      for (unsigned int gid = 0; gid < _dtlclvs.size(); ++gid) {
        _dtlclvsBackup[gid]._survivingTransferSums = _dtlclvs[gid]._survivingTransferSums;
        for (auto speciesNode: getSpeciesNodesToUpdate()) {
          auto e = speciesNode->node_index;
          _dtlclvsBackup[gid]._uq[e] = _dtlclvs[gid]._uq[e];
        }
      }
    } else { 
      std::swap(_dtlclvs, _dtlclvsBackup);
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::afterComputeLogLikelihood()
{
  AbstractReconciliationModel<REAL>::afterComputeLogLikelihood();
  if (this->_likelihoodMode == PartialLikelihoodMode::PartialSpecies) {
    if (this->_fastMode) {
      _transferExtinctionSum = _transferExtinctionSumBackup;
      for (unsigned int gid = 0; gid < _dtlclvs.size(); ++gid) {
        _dtlclvs[gid]._survivingTransferSums = _dtlclvsBackup[gid]._survivingTransferSums;
        for (auto speciesNode: getSpeciesNodesToUpdate()) {
          auto e = speciesNode->node_index;
          _dtlclvs[gid]._uq[e] = _dtlclvsBackup[gid]._uq[e];
        }
      }
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::rollbackToLastState()
{
  std::swap(_dtlclvs, _dtlclvsBackup);
}

template <class REAL>
void UndatedDTLModel<REAL>::getBestTransfer(pll_unode_t *parentGeneNode, 
  pll_rnode_t *originSpeciesNode,
  bool isVirtualRoot,
  pll_unode_t *&transferedGene,
  pll_unode_t *&stayingGene,
  pll_rnode_t *&recievingSpecies,
  REAL &proba)
{
  proba = REAL();
  auto e = originSpeciesNode->node_index;
  auto u_left = this->getLeft(parentGeneNode, isVirtualRoot);
  auto u_right = this->getRight(parentGeneNode, isVirtualRoot);
  std::unordered_set<unsigned int> parents;
  parents.insert(originSpeciesNode->node_index);
  for (auto parent = originSpeciesNode; parent->parent != 0; parent = parent->parent) {
    parents.insert(parent->parent->node_index);
  }
  for (auto species: this->_allSpeciesNodes) {
    auto h = species->node_index;
    if (parents.count(h)) {
      continue;
    }
    double factor = _PT[e] / static_cast<double>(this->_allSpeciesNodes.size());
    REAL probaLeftTransfered = (_dtlclvs[u_left->node_index]._uq[h] * _dtlclvs[u_right->node_index]._uq[e]) * factor;
    REAL probaRightTransfered = (_dtlclvs[u_right->node_index]._uq[h] * _dtlclvs[u_left->node_index]._uq[e]) * factor;
    if (proba < probaLeftTransfered) {
      proba = probaLeftTransfered;
      transferedGene = u_left;
      stayingGene = u_right;
      recievingSpecies = species;
    }
    if (proba < probaRightTransfered) {
      proba = probaRightTransfered;
      transferedGene = u_right;
      stayingGene = u_left;
      recievingSpecies = species;
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::getBestTransferLoss(Scenario &scenario,
   pll_unode_t *parentGeneNode, 
  pll_rnode_t *originSpeciesNode,
  pll_rnode_t *&recievingSpecies,
  REAL &proba)
{
  proba = REAL();
  auto e = originSpeciesNode->node_index;
  auto u = parentGeneNode->node_index;
  std::unordered_set<unsigned int> parents;
  parents.insert(originSpeciesNode->node_index);
  for (auto parent = originSpeciesNode; parent->parent != 0; parent = parent->parent) {
    parents.insert(parent->parent->node_index);
  }
  for (auto species: this->_allSpeciesNodes) {
    auto h = species->node_index;
    if (parents.count(h)) {
      continue;
    }
    REAL factor = _uE[e] * (_PT[e] / static_cast<double>(this->_allSpeciesNodes.size()));
    REAL newProba = _dtlclvs[u]._uq[h] * factor;  
    if (proba < newProba) {
      if (!scenario.isBlacklisted(u, h)) {
        scenario.blackList(u, h);
        proba = newProba;
        recievingSpecies = species;
      }
    }
  }
}
