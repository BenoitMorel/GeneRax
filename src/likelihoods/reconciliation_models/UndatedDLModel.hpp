#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>
#include <maths/ScaledValue.hpp>
#include <unordered_set>

using namespace std;

/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In this implementation, we do not allow transfers, which 
* allows a lot of algorithmic shortcuts
*/
class UndatedDLModel: public AbstractReconciliationModel {
public:
  UndatedDLModel(pll_rtree_t *speciesTree, const GeneSpeciesMapping &map);
  virtual ~UndatedDLModel();
  
  // overloaded from parent
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);  
  // overloaded from parent
  virtual void invalidateCLV(int nodeIndex);
  
protected:
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual double computeLogLikelihoodInternal(shared_ptr<pllmod_treeinfo_t> treeinfo);
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

  // ll[speciesId] = likelihood of a gene tree to be present under a species node
  vector<ScaledValue> ll; // sam
 
  // geme ids in postorder 
  vector<int> geneIds;

  // should we update all the clvs at next updateCLVs call?
  bool allCLVInvalid;

  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  unordered_set<int> invalidCLVs;

  // repeats
  vector<unsigned long> repeatsId; // repeatsId[geneId]
  vector<vector<ScaledValue> > cache_;
private:
  void getCLVsToUpdate(pllmod_treeinfo_t &treeinfo, unordered_set<int> &nodesToUpdate);
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot = false);
  void computeGeneProbabilities(pll_unode_t *geneNode, 
      vector<ScaledValue> &clv);
  void updateCLV(pll_unode_t *geneNode);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  pll_unode_t *computeLikelihoods(pllmod_treeinfo_t &treeinfo, ScaledValue &bestValue);
};

