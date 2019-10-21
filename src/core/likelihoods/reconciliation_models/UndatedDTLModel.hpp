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
  UndatedDTLModel();
  virtual ~UndatedDTLModel();
  
  // overloaded from parent
  virtual void setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates,
      const std::vector< std::vector <double > > &transferFrequencies);
protected:
  // overloaded from parent
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // overloaded from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overloaded from parent
  virtual REAL getRootLikelihood(pll_unode_t *root) const;
  // overload from parent
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot);
  virtual REAL getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {return _uq[root->node_index + this->_maxGeneId + 1][speciesRoot->node_index];}
  virtual REAL getLikelihoodFactor() const;
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false);
private:
  // model
  std::vector<double> _PD; // Duplication probability, per branch
  std::vector<double> _PL; // Loss probability, per branch
  std::vector<double> _PT; // Transfer probability, per branch
  std::vector<double> _PS; // Speciation probability, per branch

  // SPECIES
  std::vector<REAL> _uE; // Probability for a gene to become extinct on each brance
  REAL _transferExtinctionSum;
  std::vector<REAL> _ancestralExctinctionCorrection;  

  // CLVs
  // _uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  std::vector<std::vector<REAL> > _uq;
  std::vector<REAL> _survivingTransferSums;
  std::vector<std::vector<REAL> > _ancestralCorrection;



private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false);
  void updateTransferSums(REAL &transferExtinctionSum,
    std::vector<REAL> &ancestralExtinctionCorrection,
    const std::vector<REAL> &probabilities);
  void resetTransferSums(REAL &transferSum,
    std::vector<REAL> &ancestralCorrection);
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
    
  REAL getCorrectedTransferExtinctionSum(unsigned int speciesId) const {
    return (_transferExtinctionSum - _ancestralExctinctionCorrection[speciesId]) * _PT[speciesId];
  }

  REAL getCorrectedTransferSum(unsigned int geneId, unsigned int speciesId) const
  {
    return (_survivingTransferSums[geneId] - _ancestralCorrection[geneId][speciesId]) * _PT[speciesId];
  }

};


const unsigned int CACHE_SIZE = 100000;
const unsigned int IT = 5;

template <class REAL>
UndatedDTLModel<REAL>::UndatedDTLModel()
{
  this->_maxGeneId = 1;
}


template <class REAL>
void UndatedDTLModel<REAL>::setInitialGeneTree(pll_utree_t *tree)
{
  AbstractReconciliationModel<REAL>::setInitialGeneTree(tree);
  std::vector<REAL> zeros(this->speciesNodesCount_);
  _uq = std::vector<std::vector<REAL> >(2 * (this->_maxGeneId + 1),zeros);
  _survivingTransferSums = std::vector<REAL>(2 * (this->_maxGeneId + 1));
  _ancestralCorrection = std::vector<std::vector<REAL> >(2 * (this->_maxGeneId + 1),zeros);
}

template <class REAL>
void UndatedDTLModel<REAL>::updateTransferSums(REAL &transferSum,
    std::vector<REAL> &ancestralCorrection,
    const std::vector<REAL> &probabilities)
{
  transferSum = REAL();
  for (int i = static_cast<int>(this->speciesNodes_.size()) - 1; i >= 0; --i) {
    auto speciesNode = this->speciesNodes_[static_cast<unsigned int>(i)];
    auto e = speciesNode->node_index;
    ancestralCorrection[e] = probabilities[e];
    if (speciesNode->parent) {
      auto p = speciesNode->parent->node_index;
      ancestralCorrection[e] += ancestralCorrection[p];
    }
  }
  for (auto speciesNode: this->speciesNodes_) {
    auto e = speciesNode->node_index;
    ancestralCorrection[e] /= double(this->speciesNodes_.size());
    transferSum += probabilities[e];
  }
  transferSum /= this->speciesNodes_.size();
}


template <class REAL>
void UndatedDTLModel<REAL>::setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates,
      const std::vector< std::vector <double > > &)
{
  this->geneRoot_ = 0;
  assert(this->speciesNodesCount_ == dupRates.size());
  assert(this->speciesNodesCount_ == lossRates.size());
  assert(this->speciesNodesCount_ == transferRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PT = transferRates;
  _PS = std::vector<double>(this->speciesNodesCount_, 1.0);
  for (auto speciesNode: this->speciesNodes_) {
    auto e = speciesNode->node_index;
    auto sum = _PD[e] + _PL[e] + _PT[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] /= sum;
  } 
  _uE = std::vector<REAL>(this->speciesNodesCount_);
  resetTransferSums(_transferExtinctionSum, _ancestralExctinctionCorrection);
  for (unsigned int it = 0; it < IT; ++it) {
    for (auto speciesNode: this->speciesNodes_) {
      auto e = speciesNode->node_index;
      REAL proba(_PL[e]);
      proba += _uE[e] * _uE[e] * _PD[e] + getCorrectedTransferExtinctionSum(e) * _uE[e];
      if (speciesNode->left) {
        proba += _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index] * _PS[e];
      }
      ASSERT_PROBA(proba)
      _uE[speciesNode->node_index] = proba;
    }
    updateTransferSums(_transferExtinctionSum, _ancestralExctinctionCorrection, _uE);
  }
  this->invalidateAllCLVs();
}

template <class REAL>
UndatedDTLModel<REAL>::~UndatedDTLModel() { }


template <class REAL>
void UndatedDTLModel<REAL>::resetTransferSums(REAL &transferSum,
    std::vector<REAL> &ancestralCorrection)
{
  transferSum = REAL();
  if (ancestralCorrection.size()) {
    for (auto &e: ancestralCorrection) {
      e = REAL();
    }
  } else {
    ancestralCorrection = std::vector<REAL>(this->speciesNodesCount_);
  }
}


template <class REAL>
void UndatedDTLModel<REAL>::updateCLV(pll_unode_t *geneNode)
{
  auto gid = geneNode->node_index;
  for (auto speciesNode: this->speciesNodes_) {
    _uq[gid][speciesNode->node_index] = REAL();
  }
  resetTransferSums(_survivingTransferSums[gid], _ancestralCorrection[gid]);
  for (unsigned int it = 0; it < IT; ++it) {
    for (auto speciesNode: this->speciesNodes_) {
      computeProbability(geneNode, 
          speciesNode, 
          _uq[gid][speciesNode->node_index]);
    }
    updateTransferSums(_survivingTransferSums[gid], _ancestralCorrection[gid], _uq[gid]);
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
  for (auto species: this->speciesNodes_) {
    auto h = species->node_index;
    if (parents.count(h)) {
      continue;
    }
    double factor = _PT[e] / static_cast<double>(this->speciesNodes_.size());
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
  for (auto species: this->speciesNodes_) {
    auto h = species->node_index;
    if (parents.count(h)) {
      continue;
    }
    REAL factor = _uE[e] * (_PT[e] / static_cast<double>(this->speciesNodes_.size()));
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
  
  if (isSpeciesLeaf and isGeneLeaf and e == this->geneToSpecies_[gid]) {
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

  unsigned int maxValueIndex = static_cast<unsigned int>(distance(values.begin(), max_element(values.begin(), values.end())));
  // safety check
  if (values[maxValueIndex] == REAL()) {
    REAL proba;
    computeProbability(geneNode, speciesNode, proba, isVirtualRoot);
    std::cerr << "warning: null ll scenario " << _uq[gid][e] << " " << proba  << std::endl;
    assert(false);
  }
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
      scenario.addTransfer(ReconciliationEventType::EVENT_T, gid, e, transferedGene->node_index, recievingSpecies->node_index);
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


template <class REAL>
void UndatedDTLModel<REAL>::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot)
{
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  
  if (isSpeciesLeaf and isGeneLeaf and e == this->geneToSpecies_[gid]) {
    proba = REAL(_PS[e]);
    return;
  }
  
  REAL oldProba = proba;
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
  proba += oldProba * getCorrectedTransferExtinctionSum(e);
  proba += getCorrectedTransferSum(gid, e) * _uE[e];

  // DL event
  proba += oldProba * _uE[e] * (2.0 * _PD[e]); 
  //assert(proba.isProba());
}


template <class REAL>
void UndatedDTLModel<REAL>::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;;
  for (auto speciesNode: this->speciesNodes_) {
    auto e = speciesNode->node_index;
    _uq[u][e] = REAL();
  }
  resetTransferSums(_survivingTransferSums[u], _ancestralCorrection[u]);
  for (unsigned int it = 0; it < IT; ++it) {
    for (auto speciesNode: this->speciesNodes_) {
      unsigned int e = speciesNode->node_index;
      computeProbability(virtualRoot, speciesNode, _uq[u][e], true);
    }
    updateTransferSums(_survivingTransferSums[u], _ancestralCorrection[u], _uq[u]);
  }
}


template <class REAL>
REAL UndatedDTLModel<REAL>::getRootLikelihood(pll_unode_t *root) const
{
  REAL sum = REAL();
  auto u = root->node_index + this->_maxGeneId + 1;;
  for (auto speciesNode: this->speciesNodes_) {
    auto e = speciesNode->node_index;
    sum += _uq[u][e];
  }
  return sum;
}

template <class REAL>
REAL UndatedDTLModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->speciesNodes_) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - _uE[e]);
  }
  return factor;
      
}

