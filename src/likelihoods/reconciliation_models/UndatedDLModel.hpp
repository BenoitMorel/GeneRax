#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>
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
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  // overload from parent
  virtual double computeLogLikelihoodInternal(shared_ptr<pllmod_treeinfo_t> treeinfo);
  // overload from parent
  virtual void updateCLV(pll_unode_t *geneNode);
private:
  // model
  vector<double> PD; // Duplication probability, per branch
  vector<double> PL; // Loss probability, per branch
  vector<double> PS; // Speciation probability, per branch

  // SPECIES
  vector<double> uE; // Probability for a gene to become extinct on each brance
  
  // CLVs
  // uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  vector<vector<ScaledValue> > uq;



private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot = false) const;
  void computeGeneProbabilities(pll_unode_t *geneNode, 
      vector<ScaledValue> &clv);
  void computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  void updateRoot(pllmod_treeinfo_t &treeinfo);
  double getSumLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
};

