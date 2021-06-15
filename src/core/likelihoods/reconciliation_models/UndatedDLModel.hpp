#pragma once

#include <likelihoods/reconciliation_models/GTBaseReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <util/Scenario.hpp>
#include <cmath>





/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In this implementation, we do not allow transfers, which 
* allows a lot of algorithmic shortcuts
*/
template <class REAL>
class UndatedDLModel: public GTBaseReconciliationModel<REAL> {
public:
  UndatedDLModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMappingp, 
      const RecModelInfo &recModelInfo): 
    GTBaseReconciliationModel<REAL>(speciesTree, 
        geneSpeciesMappingp, 
        recModelInfo) {}
  
  
  UndatedDLModel(const UndatedDLModel &) = delete;
  UndatedDLModel & operator = (const UndatedDLModel &) = delete;
  UndatedDLModel(UndatedDLModel &&) = delete;
  UndatedDLModel & operator = (UndatedDLModel &&) = delete;
  virtual ~UndatedDLModel();
  
  // overloaded from parent
  virtual void setRates(const RatesVector &rates);
protected:
  // overload from parent
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  // overload from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overload from parent
  virtual REAL getGeneRootLikelihood(pll_unode_t *root) const;
  virtual REAL getGeneRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {
    return _dlclvs[root->node_index + this->_maxGeneId + 1][speciesRoot->node_index];
  }

  // overload from parent
  virtual void recomputeSpeciesProbabilities();
  virtual REAL getLikelihoodFactor() const;
  // overload from parent
  virtual void computeGeneRootLikelihood(pll_unode_t *virtualRoot);
  // overlead from parent
  virtual void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr,
      bool stochastic = false);
private:
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  
  // uq[cladeId][speciesId] 
  typedef std::vector<REAL> DLCLV;
  std::vector<DLCLV> _dlclvs;
 
private:
  std::vector<pll_rnode_s *> &getSpeciesNodesToUpdate() {
    //return this->_speciesNodesToUpdate;
    return this->_allSpeciesNodes;
  }

};


template <class REAL>
void UndatedDLModel<REAL>::setInitialGeneTree(PLLUnrootedTree &tree)
{
  GTBaseReconciliationModel<REAL>::setInitialGeneTree(tree);
  assert(this->_allSpeciesNodesCount);
  assert(this->_maxGeneId);
  std::vector<REAL> zeros(this->_allSpeciesNodesCount);
  _dlclvs = std::vector<std::vector<REAL> >(2 * (this->_maxGeneId + 1),zeros);
}

template <class REAL>
void UndatedDLModel<REAL>::setRates(const RatesVector &rates)
{
  assert(rates.size() == 2);
  auto &dupRates = rates[0];
  auto &lossRates = rates[1];
  assert(this->_allSpeciesNodesCount == dupRates.size());
  assert(this->_allSpeciesNodesCount == lossRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PS = std::vector<double>(this->_allSpeciesNodesCount, 1.0);
  this->_geneRoot = 0;
  for (unsigned int e = 0; e < this->_allSpeciesNodesCount; ++e) {
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  recomputeSpeciesProbabilities();
  this->invalidateAllCLVs();
  this->invalidateAllSpeciesCLVs();
}

template <class REAL>
void UndatedDLModel<REAL>::recomputeSpeciesProbabilities()
{
  if (!_uE.size()) {
    _uE = std::vector<double>(this->_allSpeciesNodesCount, 0.0);
  }
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (this->getSpeciesLeft(speciesNode)) {
      c += _PS[e] * _uE[this->getSpeciesLeft(speciesNode)->node_index]  * _uE[this->getSpeciesRight(speciesNode)->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    ASSERT_PROBA(proba)
    _uE[speciesNode->node_index] = proba;
  }
}

template <class REAL>
UndatedDLModel<REAL>::~UndatedDLModel() { }

template <class REAL>
void UndatedDLModel<REAL>::updateCLV(pll_unode_t *geneNode)
{
  assert(geneNode);
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    computeProbability(geneNode, 
        speciesNode, 
        _dlclvs[geneNode->node_index][speciesNode->node_index]);
  }
}


template <class REAL>
void UndatedDLModel<REAL>::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot,
      Scenario *,
      Scenario::Event *event,
      bool stochastic)
  
{
  auto gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto e = speciesNode->node_index;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
  }

  if (event) {
    event->geneNode = gid; 
    event->speciesNode = e;
    event->type = ReconciliationEventType::EVENT_None; 
  }


  proba = REAL();
  if (isSpeciesLeaf and isGeneLeaf) {
    // present
    if (e == this->_geneToSpecies[gid]) {
      proba = REAL(_PS[e]);
    } else {
      proba = REAL();
    }
    return;
  }
  
  typedef std::array<REAL, 8>  ValuesArray;
  ValuesArray values;
  values[0] = values[1] = values[2] = values[3] = REAL();
  values[4] = values[5] = values[6] = values[7] = REAL();
  
  if (not isGeneLeaf) {
    auto u_left = leftGeneNode->node_index;
    auto u_right = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      values[0] = _dlclvs[u_left][f];
      values[1] = _dlclvs[u_left][g];
      values[0] *= _dlclvs[u_right][g];
      values[1] *= _dlclvs[u_right][f];
      values[0] *= _PS[e]; 
      values[1] *= _PS[e]; 
      scale(values[0]);
      scale(values[1]);
      proba += values[0];
      proba += values[1];
    }
    // D event
    values[2] = _dlclvs[u_left][e];
    values[2] *= _dlclvs[u_right][e];
    values[2] *= _PD[e];
    scale(values[2]);
    proba += values[2];
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _dlclvs[gid][f];
    values[3] *= (_uE[g] * _PS[e]);
    scale(values[3]);
    values[4] = _dlclvs[gid][g];
    values[4] *=  (_uE[f] * _PS[e]);
    scale(values[4]);
    proba += values[3];
    proba += values[4];
  }
  // DL event
  proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
  //ASSERT_PROBA(proba);
  
  if (event) {
    int maxValueIndex = 0;
    if (!stochastic) {
      maxValueIndex =static_cast<int>(std::distance(values.begin(),
          std::max_element(values.begin(), values.end())
          ));
    } else {
      maxValueIndex = sampleIndex<ValuesArray, REAL>(values);
    }
    if (maxValueIndex == -1 || values[maxValueIndex] == REAL()) {
      event->type = ReconciliationEventType::EVENT_Invalid;
      return;
    }
    switch(maxValueIndex) {
    case 0:
      event->type = ReconciliationEventType::EVENT_S;
      event->leftGeneIndex = leftGeneNode->node_index;
      event->rightGeneIndex = rightGeneNode->node_index;
      break;
    case 1:
      event->type = ReconciliationEventType::EVENT_S;
      event->leftGeneIndex = rightGeneNode->node_index;
      event->rightGeneIndex = leftGeneNode->node_index;
      break;
    case 2:
      event->type = ReconciliationEventType::EVENT_D;
      event->leftGeneIndex = leftGeneNode->node_index;
      event->rightGeneIndex = rightGeneNode->node_index;
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
    default:
      assert(false);
    }
  }
}
  
template <class REAL>
REAL UndatedDLModel<REAL>::getGeneRootLikelihood(pll_unode_t *root) const
{
  REAL sum = REAL();
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    sum += _dlclvs[u][e];
  }
  return sum;
}

template <class REAL>
void UndatedDLModel<REAL>::computeGeneRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    auto e = speciesNode->node_index;
    computeProbability(virtualRoot, speciesNode, _dlclvs[u][e], true);
  }
}

template <class REAL>
REAL UndatedDLModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}


