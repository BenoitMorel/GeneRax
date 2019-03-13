#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <maths/ScaledValue.hpp>

using namespace std;




/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In this implementation, we do not allow transfers, which 
* allows a lot of algorithmic shortcuts
*/
class UndatedDLModel: public AbstractReconciliationModel {
public:
  UndatedDLModel();
  virtual ~UndatedDLModel();
  
  // overloaded from parent
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);  
protected:
  // overload from parent
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // overload from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overload from parent
  virtual ScaledValue getRootLikelihood(pll_unode_t *root) const;
  virtual ScaledValue getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {return _uq[root->node_index + _maxGeneId + 1][speciesRoot->node_index];}
  // overload from parent
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot);
  // overlead from parent
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false); //todobenoit make it pure virtual
private:
  vector<double> _PD; // Duplication probability, per species branch
  vector<double> _PL; // Loss probability, per species branch
  vector<double> _PS; // Speciation probability, per species branch
  vector<double> _uE; // Extinction probability, per species branch
  
  // uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  vector<vector<ScaledValue> > _uq;

private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot = false) const;
};

