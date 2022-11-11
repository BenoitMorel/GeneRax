#pragma once

#include <likelihoods/reconciliation_models/GTBaseReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <util/Scenario.hpp>
#include <cmath>




/**
 * TOOD DOCUMENTATION: check what this model does, and where it is used
 */
template <class REAL>
class SimpleDSModel: public GTBaseReconciliationModel<REAL> {
public:
  SimpleDSModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMappingp, 
      const RecModelInfo &recModelInfo): 
    GTBaseReconciliationModel<REAL>(speciesTree, 
        geneSpeciesMappingp, 
        recModelInfo) {}
  
  
  SimpleDSModel(const SimpleDSModel &) = delete;
  SimpleDSModel & operator = (const SimpleDSModel &) = delete;
  SimpleDSModel(SimpleDSModel &&) = delete;
  SimpleDSModel & operator = (SimpleDSModel &&) = delete;
  virtual ~SimpleDSModel();
  
  // overloaded from parent
  virtual void setRates(const RatesVector &rates);
protected:
  // overload from parent
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  // overload from parent
  virtual void updateCLV(corax_unode_t *geneNode);
  // overload from parent
  virtual REAL getGeneRootLikelihood(corax_unode_t *root) const;
  virtual REAL getGeneRootLikelihood(corax_unode_t *root, corax_rnode_t *) {
    return _dsclvs[root->node_index].proba;
  }

  // overload from parent
  virtual REAL getLikelihoodFactor() const {return REAL(1.0);}
  // overload from parent
  virtual void recomputeSpeciesProbabilities(){};
  // overload from parent
  virtual void computeGeneRootLikelihood(corax_unode_t *virtualRoot);
  // overlead from parent
  virtual void computeProbability(corax_unode_t *geneNode, corax_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr,
      bool stochastic = false);
private:
  double _PS; // Speciation probability
  double _PD; // Duplication probability
  
  struct DSCLV {
    REAL proba;
    std::set<unsigned int> clade;
    unsigned int genesCount;
    DSCLV():proba(REAL()), genesCount(0) {}
  };
  std::vector<DSCLV> _dsclvs;
};


template <class REAL>
void SimpleDSModel<REAL>::setInitialGeneTree(PLLUnrootedTree &tree)
{
  GTBaseReconciliationModel<REAL>::setInitialGeneTree(tree);
  assert(this->_maxGeneId);
  _dsclvs = std::vector<DSCLV>(2 * (this->_maxGeneId + 1));
}

template <class REAL>
void SimpleDSModel<REAL>::setRates(const RatesVector &rates)
{
  assert(rates.size() == 1);
  auto &dupRates = rates[0];
  _PD = dupRates[0];
  _PS = 1.0;
  auto sum = _PD + _PS;
  _PD /= sum;
  _PS /= sum;
  this->_geneRoot = 0;
  this->invalidateAllCLVs();
}


template <class REAL>
SimpleDSModel<REAL>::~SimpleDSModel() { }

template <class REAL>
void SimpleDSModel<REAL>::updateCLV(corax_unode_t *geneNode)
{
  assert(geneNode);
  computeProbability(geneNode, 
      nullptr, 
      _dsclvs[geneNode->node_index].proba);
}

template <class REAL>
static REAL dividePowerTwo(REAL v, unsigned int powerTwo)
{
  /*
  v *= pow(2.0, -double(powerTwo));
  return v;
  */
  const unsigned int MAX_EXPO_TWO = 16;
  const double MAX_POWER_TWO = pow(2.0, -double(MAX_EXPO_TWO));
  
  while (powerTwo > MAX_EXPO_TWO) {
    v *= MAX_POWER_TWO;
    scale<REAL>(v);
    powerTwo -= MAX_EXPO_TWO;
  }
  v *= pow(2.0, -double(powerTwo));
  scale<REAL>(v);
  return v;
}

template <class REAL>
void SimpleDSModel<REAL>::computeProbability(corax_unode_t *geneNode, 
    corax_rnode_t *, 
      REAL &proba,
      bool isVirtualRoot,
      Scenario *,
      Scenario::Event *event,
      bool)
  
{
  assert(!event); // we cannot reconcile with this model
  auto gid = geneNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  auto &clade = _dsclvs[gid].clade;

  if (isGeneLeaf) {
    auto speciesId = this->_geneToSpecies[gid];
    proba = REAL(_PS);
    clade.insert(speciesId);
    _dsclvs[gid].genesCount = 1;
    return;
  }
  corax_unode_t *leftGeneNode = 0;     
  corax_unode_t *rightGeneNode = 0;     
  leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
  rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  auto v = leftGeneNode->node_index;
  auto w = rightGeneNode->node_index;
  auto &leftClade = _dsclvs[v].clade;
  auto &rightClade = _dsclvs[w].clade;
  clade.clear();
  clade.insert(leftClade.begin(), leftClade.end());
  clade.insert(rightClade.begin(), rightClade.end());
  _dsclvs[gid].genesCount = _dsclvs[v].genesCount + _dsclvs[w].genesCount;
  proba = REAL(_dsclvs[v].proba * _dsclvs[w].proba);
  if (clade.size() == leftClade.size() + rightClade.size()) {
    // proba *=  _PS / 2^(species - 1)
    proba *= _PS;
    unsigned int species = clade.size();
    proba = dividePowerTwo(proba, species - 1);
  } else {
    // proba *= _PD * (2^(genes - 1) - 2^(species -1))
    // or proba *= _PD / ((2^(species - 1)) * (2^(genes - species) - 1))
    unsigned int species = clade.size();
    unsigned int genes = _dsclvs[gid].genesCount;
    proba *= _PD;
    proba = dividePowerTwo(proba, species - 1);
    unsigned int diff = genes - species;
    if (diff < 8) { // no overflow, and we account for "- 1.0"
      proba /= (pow(2.0, diff) - 1.0);
      scale<REAL>(proba);
    } else { // we neglect the "- 1.0" and make sure we do not overflow
      proba = dividePowerTwo(proba, diff);
    }
  }
    //proba /= pow(2.0, _dsclvs[gid].genesCount - 1) - pow(2.0, clade.size() - 1);
}
  
template <class REAL>
REAL SimpleDSModel<REAL>::getGeneRootLikelihood(corax_unode_t *root) const
{
  auto u = root->node_index + this->_maxGeneId + 1;
  return _dsclvs[u].proba;
}

template <class REAL>
void SimpleDSModel<REAL>::computeGeneRootLikelihood(corax_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  computeProbability(virtualRoot, nullptr, _dsclvs[u].proba, true);
}



