#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>

static const int TRANSFER_FREQUENCIES_EPSILON = 0.1;

/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In addition, we forbid transfers to parent species
*/
template <class REAL>
class UndatedDTLModelAdvanced: public AbstractReconciliationModel<REAL> {
public:
  UndatedDTLModelAdvanced();
  virtual ~UndatedDTLModelAdvanced();
  
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
  std::vector<std::vector<double> > _PT; // Transfer probability, per pair of branches
  std::vector<double> _PS; // Speciation probability, per branch

  // SPECIES
  std::vector<REAL> _uE; // Probability for a gene to become extinct on each brance

  // CLVs
  // _uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  std::vector<std::vector<REAL> > _uq;



private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false);
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

};


const unsigned int IT_ADVANCED = 5;

template <class REAL>
UndatedDTLModelAdvanced<REAL>::UndatedDTLModelAdvanced()
{
  this->_maxGeneId = 1;
}


template <class REAL>
void UndatedDTLModelAdvanced<REAL>::setInitialGeneTree(pll_utree_t *tree)
{
  AbstractReconciliationModel<REAL>::setInitialGeneTree(tree);
  std::vector<REAL> zeros(this->speciesNodesCount_);
  _uq = std::vector<std::vector<REAL> >(2 * (this->_maxGeneId + 1),zeros);
}

template <class REAL>
void UndatedDTLModelAdvanced<REAL>::setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates,
      const std::vector< std::vector <double > > &transferFrequencies)
{
  this->geneRoot_ = 0;
  assert(this->speciesNodesCount_ == dupRates.size());
  assert(this->speciesNodesCount_ == lossRates.size());
  assert(this->speciesNodesCount_ == transferRates.size());
  assert(this->speciesNodesCount_ == transferFrequencies.size());
  _PD = dupRates;
  _PL = lossRates;
  _PT = std::vector<std::vector<double> >(this->speciesNodesCount_, transferRates);
  _PS = std::vector<double>(this->speciesNodesCount_, 1.0);
  auto normalizedTransferFrequencies = transferFrequencies;
  for (unsigned int e = 0; e < this->speciesNodesCount_; ++e) {
    double sum = 0.0;
    assert(normalizedTransferFrequencies[e].size() == this->speciesNodesCount_);
    for (unsigned int f = 0; f < this->speciesNodesCount_; ++f) {
      normalizedTransferFrequencies[e][f] += TRANSFER_FREQUENCIES_EPSILON;
      sum += normalizedTransferFrequencies[e][f];
    }
    for (unsigned int f = 0; f < this->speciesNodesCount_; ++f) {
      normalizedTransferFrequencies[e][f] /= sum;
    }
  }
  for (unsigned int e = 0; e < this->speciesNodesCount_; ++e) {
    auto sum = _PD[e] + _PL[e] + _PS[e];
    for (unsigned int f = 0; f < this->speciesNodesCount_; ++f) {
      _PT[e][f] *= normalizedTransferFrequencies[e][f];
      sum += _PT[e][f];
    }
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
    for (auto &PT: _PT[e]) {
      PT /= sum;
    }
  } 
  _uE = std::vector<REAL>(this->speciesNodesCount_);
  for (unsigned int it = 0; it < IT_ADVANCED; ++it) {
    for (auto speciesNode: this->speciesNodes_) {
      auto e = speciesNode->node_index;
      REAL proba(_PL[e]);
      proba += _uE[e] * _uE[e] * _PD[e];// + getTransferExtinctionSum(e) * _uE[e];
      auto transferTerm = REAL();
      for (auto recievingSpeciesNode: this->speciesNodes_) {
        auto h = recievingSpeciesNode->node_index;
        transferTerm += _uE[h] * _PT[e][h];
      }
      proba += transferTerm * _uE[e] * (1 / double(this->speciesNodesCount_));
      if (speciesNode->left) {
        proba += _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index] * _PS[e];
      }
      ASSERT_PROBA(proba)
      _uE[speciesNode->node_index] = proba;
    }
  }
  this->invalidateAllCLVs();
}

template <class REAL>
UndatedDTLModelAdvanced<REAL>::~UndatedDTLModelAdvanced() { }




template <class REAL>
void UndatedDTLModelAdvanced<REAL>::updateCLV(pll_unode_t *geneNode)
{
  auto gid = geneNode->node_index;
  for (auto speciesNode: this->speciesNodes_) {
    _uq[gid][speciesNode->node_index] = REAL();
  }
  for (unsigned int it = 0; it < IT_ADVANCED; ++it) {
    for (auto speciesNode: this->speciesNodes_) {
      computeProbability(geneNode, 
          speciesNode, 
          _uq[gid][speciesNode->node_index]);
    }
  }
}



template <class REAL>
void UndatedDTLModelAdvanced<REAL>::backtrace(pll_unode_t *, pll_rnode_t *, 
      Scenario &,
      bool ) 
{
  Logger::info << "BACKTRACE NOT IMPLEMENTED" << std::endl;
}


template <class REAL>
void UndatedDTLModelAdvanced<REAL>::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
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
    REAL transferTerm = REAL();
    for (auto *recievingSpeciesNode: this->speciesNodes_) {
      auto h = recievingSpeciesNode->node_index;
      transferTerm += (_uq[u_left][h] * _uq[u_right][e] + _uq[u_left][e] * _uq[u_right][h]) * (_PT[e][h] / this->speciesNodesCount_) ;
    }
    proba += transferTerm;
  }
  if (not isSpeciesLeaf) {
    // SL event
    proba += (_uq[gid][f] * _uE[g] + _uq[gid][g] * _uE[f]) *_PS[e];
  }
  // TL event
  REAL transferLossTerm = REAL();
  for (auto *recievingSpeciesNode: this->speciesNodes_) {
    auto h = recievingSpeciesNode->node_index;
    transferLossTerm += (oldProba * _uE[h] + _uE[e] * _uq[gid][h]) * _PT[e][h];
  }
  proba += transferLossTerm;

  // DL event
  proba += oldProba * _uE[e] * (2.0 * _PD[e]); 
}


template <class REAL>
void UndatedDTLModelAdvanced<REAL>::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;;
  for (auto speciesNode: this->speciesNodes_) {
    auto e = speciesNode->node_index;
    _uq[u][e] = REAL();
  }
  for (unsigned int it = 0; it < IT_ADVANCED; ++it) {
    for (auto speciesNode: this->speciesNodes_) {
      unsigned int e = speciesNode->node_index;
      computeProbability(virtualRoot, speciesNode, _uq[u][e], true);
    }
  }
}


template <class REAL>
REAL UndatedDTLModelAdvanced<REAL>::getRootLikelihood(pll_unode_t *root) const
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
REAL UndatedDTLModelAdvanced<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->speciesNodes_) {
    auto e = speciesNode->node_index;
    if (!(_uE[e] < REAL(1.0))) {
      std::cerr << "error : " <<  _uE[e] << std::endl;
    }
    factor += (REAL(1.0) - _uE[e]);
  }
  return factor;
      
}

