#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>

#define IS_PROBA(x) (REAL(0.0) <= (x) && (x) <= REAL(1.0))
  
#define PRINT_ERROR_PROBA(x)  if (!IS_PROBA(proba)) {std::cerr << "error " << proba << std::endl;} assert(IS_PROBA(x));  
#define ONEMORE
#define MYINVARIANT

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
protected:
  // overloaded from parent
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // overloaded from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overloaded from parent
  virtual REAL getRootLikelihood(pll_unode_t *root) const;
  // overload from parent
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot);
  virtual REAL getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {
    return _uq[root->node_index + this->_maxGeneId + 1][speciesRoot->node_index];
  }
  virtual REAL getLikelihoodFactor() const;
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false);
  virtual void beforeComputeLogLikelihood(); 
  virtual void afterComputeLogLikelihood(); 
private:
  // model
  std::vector<double> _PD; // Duplication probability, per branch
  std::vector<double> _PL; // Loss probability, per branch
  std::vector<double> _PT; // Transfer probability, per branch
  std::vector<double> _PS; // Speciation probability, per branch
  // SPECIES
  std::vector<REAL> _uE; // Probability for a gene to become extinct on each brance
  REAL _transferExtinctionSum;

  // CLVs
  // _uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  std::vector<std::vector<REAL> > _uq;
  std::vector<std::vector<REAL> > _uqBackup;
  std::vector<REAL> _survivingTransferSums;
  
  // fast mode 
  REAL _transferExtinctionSumBackup;
  std::vector<REAL> _survivingTransferSumsBackup;
  std::vector<REAL> _survivingTransferSumsInvariant;
  std::vector<REAL> _survivingTransferSumsOneMore;

private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false);
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
  void getBestTransferLoss(pll_unode_t *parentGeneNode, 
    pll_rnode_t *originSpeciesNode,
    pll_rnode_t *&recievingSpecies,
    REAL &proba);
  unsigned int getIterationsNumber() const { return this->_fastMode ? 1 : 5;}    
  REAL getCorrectedTransferExtinctionSum(unsigned int speciesId) const {
    return _transferExtinctionSum * _PT[speciesId];
  }

  REAL getCorrectedTransferSum(unsigned int geneId, unsigned int speciesId) const
  {
    return _survivingTransferSums[geneId] * _PT[speciesId];
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
  _uq = std::vector<std::vector<REAL> >(2 * (this->_maxGeneId + 1),zeros);
  _uqBackup = std::vector<std::vector<REAL> >(2 * (this->_maxGeneId + 1),zeros);
  _survivingTransferSums = std::vector<REAL>(2 * (this->_maxGeneId + 1));
  _survivingTransferSumsOneMore = std::vector<REAL>(2 * (this->_maxGeneId + 1));
  _survivingTransferSumsBackup = std::vector<REAL>(2 * (this->_maxGeneId + 1));
  _survivingTransferSumsInvariant = std::vector<REAL>(2 * (this->_maxGeneId + 1));
}

  template <class REAL>
void UndatedDTLModel<REAL>::resetTransferSums(const REAL &transferSum,
    REAL &transferSumInvariant,
    const std::vector<REAL> &probabilities)
{
#ifdef MYINVARIANT
  if (this->_fastMode) {
    REAL diff = REAL();
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      diff += probabilities[speciesNode->node_index];
    }
    diff /= this->_allSpeciesNodes.size();
    transferSumInvariant = transferSum - diff;
  }
#endif
}

template <class REAL>
void UndatedDTLModel<REAL>::updateTransferSums(REAL &transferSum,
    const REAL &transferSumInvariant,
    const std::vector<REAL> &probabilities)
{
  transferSum = REAL();
#ifdef MYINVARIANT
  for (auto speciesNode:  getSpeciesNodesToUpdate()) {
#else 
  for (auto speciesNode: this->_allSpeciesNodes) {
#endif
    auto e = speciesNode->node_index;
    transferSum += probabilities[e];
  }
  transferSum /= this->_allSpeciesNodes.size();
#ifdef MYINVARIANT 
  if (this->_fastMode) {
    transferSum += transferSumInvariant;
  }
#endif
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
  _PS = std::vector<double>(this->_allSpeciesNodesCount, 1.0);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    auto sum = _PD[e] + _PL[e] + _PT[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] /= sum;
  } 
  _uE = std::vector<REAL>(this->_allSpeciesNodesCount);
  REAL unused = REAL(1.0);
  resetTransferSums(_transferExtinctionSum, unused, _uE);
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    for (auto speciesNode: this->_allSpeciesNodes) {
      auto e = speciesNode->node_index;
      REAL proba(_PL[e]);
      proba += _uE[e] * _uE[e] * _PD[e] + getCorrectedTransferExtinctionSum(e) * _uE[e];
      if (speciesNode->left) {
        proba += _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index] * _PS[e];
      }
      //PRINT_ERROR_PROBA(proba)
      _uE[speciesNode->node_index] = proba;
    }
    updateTransferSums(_transferExtinctionSum, unused, _uE);
  }
  this->invalidateAllCLVs();
  this->invalidateAllSpeciesCLVs();
}

template <class REAL>
UndatedDTLModel<REAL>::~UndatedDTLModel() { }




template <class REAL>
void UndatedDTLModel<REAL>::updateCLV(pll_unode_t *geneNode)
{
  auto gid = geneNode->node_index;
#ifdef ONEMORE
  resetTransferSums(this->_fastMode ? _survivingTransferSumsOneMore[gid] : _survivingTransferSums[gid], _survivingTransferSumsInvariant[gid], _uq[gid]);
#else
  resetTransferSums(_survivingTransferSums[gid], _survivingTransferSumsInvariant[gid], _uq[gid]);
#endif
  
  if (!this->_fastMode) {
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      _uq[gid][speciesNode->node_index] = REAL();
    }
  }
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    updateTransferSums(_survivingTransferSums[gid], _survivingTransferSumsInvariant[gid], _uq[gid]);
    for (auto speciesNode: getSpeciesNodesToUpdate()) { 
      computeProbability(geneNode, 
          speciesNode, 
          _uq[gid][speciesNode->node_index]);
    }
  }
  if (!this->_fastMode) {
    updateTransferSums(_survivingTransferSumsOneMore[gid], _survivingTransferSumsInvariant[gid], _uq[gid]);
  }
}



template <class REAL>
void UndatedDTLModel<REAL>::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot)
{
  
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  
  if (isSpeciesLeaf and isGeneLeaf and e == this->_geneToSpecies[gid]) {
    proba = REAL(_PS[e]);
    return;
  }
  
  //REAL oldProba = proba;
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
      proba += (_uq[u_left][f] * _uq[u_right][g] + _uq[u_left][g] * _uq[u_right][f]) * _PS[e];
    }
    // D event
    REAL temp = _uq[u_left][e];
    temp *= _uq[u_right][e];
    temp *= _PD[e];
    proba += temp;
    // T event
    proba += getCorrectedTransferSum(u_left, e) * _uq[u_right][e]; 
    proba += getCorrectedTransferSum(u_right, e) * _uq[u_left][e]; 
  }
  if (not isSpeciesLeaf) {
    // SL event
    proba += (_uq[gid][f] * _uE[g] + _uq[gid][g] * _uE[f]) *_PS[e];
  }
  // TL event
  //proba += oldProba * getCorrectedTransferExtinctionSum(e);
  proba += getCorrectedTransferSum(gid, e) * _uE[e];

  // DL event
  //proba += oldProba * _uE[e] * (2.0 * _PD[e]); 
  //assert(IS_PROBA(proba));
}


template <class REAL>
void UndatedDTLModel<REAL>::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
#ifdef ONEMORE
  resetTransferSums(this->_fastMode ? _survivingTransferSumsOneMore[u] : _survivingTransferSums[u], _survivingTransferSumsInvariant[u], _uq[u]);
#else
  resetTransferSums(_survivingTransferSums[u], _survivingTransferSumsInvariant[u], _uq[u]);
#endif
  if (!this->_fastMode) {
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      auto e = speciesNode->node_index;
      _uq[u][e] = REAL();
    }
  }
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    updateTransferSums(_survivingTransferSums[u], _survivingTransferSumsInvariant[u], _uq[u]);
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      unsigned int e = speciesNode->node_index;
      computeProbability(virtualRoot, speciesNode, _uq[u][e], true);
    }
  }
  if (!this->_fastMode) {
    updateTransferSums(_survivingTransferSumsOneMore[u], _survivingTransferSumsInvariant[u], _uq[u]);
  }
}


template <class REAL>
REAL UndatedDTLModel<REAL>::getRootLikelihood(pll_unode_t *root) const
{
  REAL sum = REAL();
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    sum += _uq[u][e];
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
void UndatedDTLModel<REAL>::beforeComputeLogLikelihood()
{
  AbstractReconciliationModel<REAL>::beforeComputeLogLikelihood();
  if (this->_fastMode) {
    _transferExtinctionSumBackup = _transferExtinctionSum;
    _survivingTransferSumsBackup = _survivingTransferSums;
    for (unsigned int gid = 0; gid < _uq.size(); ++gid) {
      for (auto speciesNode: getSpeciesNodesToUpdate()) {
        auto e = speciesNode->node_index;
        _uqBackup[gid][e] = _uq[gid][e];
      }
    }
  }
}

  template <class REAL>
void UndatedDTLModel<REAL>::afterComputeLogLikelihood()
{
  AbstractReconciliationModel<REAL>::afterComputeLogLikelihood();
  if (this->_fastMode) {
    _transferExtinctionSum = _transferExtinctionSumBackup;
    _survivingTransferSums = _survivingTransferSumsBackup;
    for (unsigned int gid = 0; gid < _uq.size(); ++gid) {
      for (auto speciesNode: getSpeciesNodesToUpdate()) {
        auto e = speciesNode->node_index;
        _uq[gid][e] = _uqBackup[gid][e];
      }
    }
  }
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
    REAL probaLeftTransfered = (_uq[u_left->node_index][h] * _uq[u_right->node_index][e]) * factor;
    REAL probaRightTransfered = (_uq[u_right->node_index][h] * _uq[u_left->node_index][e]) * factor;
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
void UndatedDTLModel<REAL>::getBestTransferLoss(pll_unode_t *parentGeneNode, 
  pll_rnode_t *originSpeciesNode,
  pll_rnode_t *&recievingSpecies,
  REAL &proba)
{
  proba = REAL();
  auto e = originSpeciesNode->node_index;
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
    REAL newProba = _uq[parentGeneNode->node_index][h] * factor;  
    if (proba < newProba) {
      proba = newProba;
      recievingSpecies = species;
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot) 
{
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  
  if (isSpeciesLeaf and isGeneLeaf and e == this->_geneToSpecies[gid]) {
    scenario.addEvent(ReconciliationEventType::EVENT_None, gid, e);
    return;
  }
  
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }
  unsigned int f = 0;
  unsigned int g = 0;
  unsigned int u_left = 0;
  unsigned int u_right = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  
  std::vector<REAL> values(8);
  // for transfers only:
  pll_unode_t *transferedGene = 0;
  pll_unode_t *stayingGene = 0;
  pll_rnode_t *recievingSpecies = 0; 
  pll_rnode_t *tlRecievingSpecies = 0; 
  if (not isGeneLeaf) {
    // S event
    u_left = leftGeneNode->node_index;
    u_right = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      values[0] = _uq[u_left][f] * _uq[u_right][g] * _PS[e];
      values[1] = _uq[u_left][g] * _uq[u_right][f] * _PS[e];
    }
    // D event
    values[2] = _uq[u_left][e];
    values[2] *= _uq[u_right][e];
    values[2] *= _PD[e];
    // T event
    getBestTransfer(geneNode, speciesNode, isVirtualRoot, transferedGene, stayingGene, recievingSpecies, values[5]);
    //values[5] = getCorrectedTransferSum(u_left, e) * _uq[u_right][e]; 
    //values[6] = getCorrectedTransferSum(u_right, e) * _uq[u_left][e]; 
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _uq[gid][f] * _uE[g] * _PS[e];
    values[4] = _uq[gid][g] * _uE[f] * _PS[e];
  }
  getBestTransferLoss(geneNode, speciesNode, tlRecievingSpecies, values[6]);

  unsigned int maxValueIndex = static_cast<unsigned int>(distance(values.begin(), 
        max_element(values.begin(), values.end())
        ));
  // safety check
  assert(!(values[maxValueIndex] <= REAL()));
  switch(maxValueIndex) {
    case 0: 
      scenario.addEvent(ReconciliationEventType::EVENT_S, gid, e);
      backtrace(leftGeneNode, speciesNode->left, scenario); 
      backtrace(rightGeneNode, speciesNode->right, scenario); 
      break;
    case 1:
      scenario.addEvent(ReconciliationEventType::EVENT_S, gid, e);
      backtrace(leftGeneNode, speciesNode->right, scenario); 
      backtrace(rightGeneNode, speciesNode->left, scenario); 
      break;
    case 2:
      scenario.addEvent(ReconciliationEventType::EVENT_D, gid, e);
      backtrace(leftGeneNode, speciesNode, scenario); 
      backtrace(rightGeneNode, speciesNode, scenario); 
      break;
    case 3: 
      scenario.addEvent(ReconciliationEventType::EVENT_SL, gid, e, speciesNode->left->node_index);
      backtrace(geneNode, speciesNode->left, scenario); 
      break;
    case 4:
      scenario.addEvent(ReconciliationEventType::EVENT_SL, gid, e, speciesNode->right->node_index);
      backtrace(geneNode, speciesNode->right, scenario); 
      break;
    case 5:
      assert(transferedGene);
      assert(recievingSpecies);
      scenario.addTransfer(ReconciliationEventType::EVENT_T, gid, e, 
          transferedGene->node_index, 
          recievingSpecies->node_index);
      backtrace(transferedGene, recievingSpecies, scenario);
      backtrace(stayingGene, speciesNode, scenario);
      break;
    case 6:
      assert(tlRecievingSpecies);
      scenario.addTransfer(ReconciliationEventType::EVENT_TL, gid, e, gid, tlRecievingSpecies->node_index); 
      backtrace(geneNode, tlRecievingSpecies, scenario); 
      break;
    default:
      std::cerr << "event " << maxValueIndex << std::endl;
      Logger::error << "Invalid event in UndatedDLModel::backtrace" << std::endl;
      assert(false);
      break;
  }

}


