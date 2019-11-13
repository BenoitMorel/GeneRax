#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
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
class UndatedDLModel: public AbstractReconciliationModel<REAL> {
public:
  UndatedDLModel(PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMappingp, bool rootedGeneTree):
    AbstractReconciliationModel<REAL>(speciesTree, geneSpeciesMappingp, rootedGeneTree) {}
  
  UndatedDLModel(const UndatedDLModel &) = delete;
  UndatedDLModel & operator = (const UndatedDLModel &) = delete;
  UndatedDLModel(UndatedDLModel &&) = delete;
  UndatedDLModel & operator = (UndatedDLModel &&) = delete;
  virtual ~UndatedDLModel();
  
  // overloaded from parent
  virtual void setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates);
protected:
  // overload from parent
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // overload from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overload from parent
  virtual REAL getRootLikelihood(pll_unode_t *root) const;
  virtual REAL getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {
    return _uq[root->node_index + this->_maxGeneId + 1][speciesRoot->node_index];
  }
  virtual REAL getLikelihoodFactor() const;
  // overload from parent
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot);
  // overlead from parent
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false); //todobenoit make it pure virtual
private:
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  
  // uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  std::vector<std::vector<REAL> > _uq;

private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false) const;
  void accountForSpeciesRoot(pll_unode_t *virtualRoot);
};


template <class REAL>
void UndatedDLModel<REAL>::setInitialGeneTree(pll_utree_t *tree)
{
  AbstractReconciliationModel<REAL>::setInitialGeneTree(tree);
  assert(this->_speciesNodesCount);
  assert(this->_maxGeneId);
  std::vector<REAL> zeros(this->_speciesNodesCount);
  _uq = std::vector<std::vector<REAL> >(2 * (this->_maxGeneId + 1),zeros);
}

static double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

template <class REAL>
void UndatedDLModel<REAL>::setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates)
{
  assert(this->_speciesNodesCount == dupRates.size());
  assert(this->_speciesNodesCount == lossRates.size());
  (void)transferRates;  
  _PD = dupRates;
  _PL = lossRates;
  _PS = std::vector<double>(this->_speciesNodesCount, 1.0);
  this->_geneRoot = 0;
  for (auto speciesNode: this->_speciesNodes) {
    auto e = speciesNode->node_index;
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  _uE = std::vector<double>(this->_speciesNodesCount, 0.0);
  for (auto speciesNode: this->_speciesNodes) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (speciesNode->left) {
      c += _PS[e] * _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    ASSERT_PROBA(proba)
    _uE[speciesNode->node_index] = proba;
  }
  this->invalidateAllCLVs();
}

template <class REAL>
UndatedDLModel<REAL>::~UndatedDLModel() { }

template <class REAL>
void UndatedDLModel<REAL>::updateCLV(pll_unode_t *geneNode)
{
  assert(geneNode);
  for (auto speciesNode: this->_speciesNodes) {
    computeProbability(geneNode, 
        speciesNode, 
        _uq[geneNode->node_index][speciesNode->node_index]);
  }
}


template <class REAL>
void UndatedDLModel<REAL>::backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot) 
{
  assert(geneNode);
  assert(speciesNode);
  auto gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;   
  std::vector<REAL> values(5, REAL());
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !speciesNode->left;
  auto e = speciesNode->node_index;
  unsigned int f = 0;
  unsigned int g = 0;
  assert(!(_uq[gid][e] <= REAL())); // check that this scenario is possible
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  if (isSpeciesLeaf and isGeneLeaf and e == this->_geneToSpecies[gid]) {
    // present
    scenario.addEvent(ReconciliationEventType::EVENT_None, gid, e);
    return;
  }
  if (not isGeneLeaf) {
    auto gp_i = leftGeneNode->node_index;
    auto gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      values[0] = _uq[gp_i][f] * _uq[gpp_i][g] * _PS[e];
      values[1] = _uq[gp_i][g] * _uq[gpp_i][f] * _PS[e];
    }
    // D event
    values[2] = _uq[gp_i][e];
    values[2] *= _uq[gpp_i][e];
    values[2] *= _PD[e];
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _uq[gid][f] * (_uE[g] * _PS[e]);
    values[4] = _uq[gid][g] * (_uE[f] * _PS[e]);
  }

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
    default:
      std::cerr << "event " << maxValueIndex << std::endl;
      Logger::error << "Invalid event in UndatedDLModel<REAL>::backtrace" << std::endl;
      assert(false);
      break;
  }
}

template <class REAL>
void UndatedDLModel<REAL>::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot) const
{
  auto gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !speciesNode->left;
  auto e = speciesNode->node_index;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
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
  if (not isGeneLeaf) {
    auto gp_i = leftGeneNode->node_index;
    auto gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      proba += (_uq[gp_i][f] * _uq[gpp_i][g] + _uq[gp_i][g] * _uq[gpp_i][f]) *_PS[e];
    }
    // D event
    REAL temp = _uq[gp_i][e];
    temp *= _uq[gpp_i][e];
    temp *= _PD[e];
    proba += temp;
  }
  if (not isSpeciesLeaf) {
    // SL event
    proba += (_uq[gid][f] * _uE[g] + _uq[gid][g] * _uE[f]) * _PS[e];
  }
  // DL event
  proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
  ASSERT_PROBA(proba);
}
  
template <class REAL>
REAL UndatedDLModel<REAL>::getRootLikelihood(pll_unode_t *root) const
{
  REAL sum = REAL();
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_speciesNodes) {
    auto e = speciesNode->node_index;
    sum += _uq[u][e];
  }
  return sum;
}

template <class REAL>
void UndatedDLModel<REAL>::accountForSpeciesRoot(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  std::vector<REAL> save_uq(_uq[u]);
  auto leftGeneNode = this->getLeft(virtualRoot, true);
  auto rightGeneNode = this->getRight(virtualRoot, true);
  for (auto speciesNode: this->_speciesNodes) {
    auto e = speciesNode->node_index;
    REAL proba;
    // D
      auto gp_i = leftGeneNode->node_index;
      auto gpp_i = rightGeneNode->node_index;
      REAL temp = _uq[gp_i][e];
      temp *= _uq[gpp_i][e];
      temp *= _PD[e];
      proba += temp;
    // no event
    proba += save_uq[e] * (1.0 - _PD[e]);
    // DL
    proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
    
    _uq[u][e] = proba;
  }
}

template <class REAL>
void UndatedDLModel<REAL>::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  for (auto speciesNode: this->_speciesNodes) {
    auto e = speciesNode->node_index;
    computeProbability(virtualRoot, speciesNode, _uq[u][e], true);
  }
//  accountForSpeciesRoot(virtualRoot);
}

template <class REAL>
REAL UndatedDLModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->_speciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}


