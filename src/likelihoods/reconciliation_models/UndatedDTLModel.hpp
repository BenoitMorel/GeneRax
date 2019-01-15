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
class UndatedDTLModel: public AbstractReconciliationModel {
public:
  UndatedDTLModel();
  virtual ~UndatedDTLModel();
  
  // overloaded from parent
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);  
  // overloaded from parent
  virtual void invalidateCLV(int nodeIndex);
    
  virtual bool implementsTransfers() {return true;}
    
protected:
  virtual void setSpeciesTree(pll_rtree_t *speciesTree);
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual double computeLogLikelihoodInternal(shared_ptr<pllmod_treeinfo_t> treeinfo);
private:
  // model
  vector<double> PD; // Duplication probability, per branch
  vector<double> PL; // Loss probability, per branch
  vector<double> PT; // Transfer probability, per branch
  vector<double> PS; // Speciation probability, per branch

  // SPECIES
  vector<double> uE; // Probability for a gene to become extinct on each brance
  double transferExtinctionSum;
  vector<vector<pll_rnode_t *> > invalidTransfers; // todobenoit this might not be the optimal implementation

  // CLVs
  // uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  vector<vector<ScaledValue> > uq;
  vector<ScaledValue> survivingTransferSums;
  vector<vector<ScaledValue> > ancestralCorrection;

  // ll[speciesId] = likelihood of the (rooted or unrooted) gene tree to be present under a species node
  vector<ScaledValue> ll; // sam
 
  // geme ids in postorder 
  vector<int> geneIds;

  // should we update all the clvs at next updateCLVs call?
  bool allCLVInvalid;

  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  unordered_set<int> invalidatedNodes;

  // is the CLV up to date?
  vector<bool> isCLVUpdated;

  int maxId;

  // repeats
  vector<unsigned long> repeatsId; // repeatsId[geneId]
  vector<vector<ScaledValue> > cache_;
private:
  void getCLVsToUpdate(pllmod_treeinfo_t &treeinfo);
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot = false) const;
  void computeGeneProbabilities(pll_unode_t *geneNode);
  void updateCLV(pll_unode_t *geneNode);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  void computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  void updateRoot(pllmod_treeinfo_t &treeinfo);
  double getSumLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  void invalidateAllCLVs();
  void updateCLVsRec(pll_unode_t *node);
  void markInvalidatedNodes(pllmod_treeinfo_t &treeinfo);
  void markInvalidatedNodesRec(pll_unode_t *node);
  void updateExtinctionTransferSums(vector<double> &ancestralExctinctionCorrection);
  void updateTransferSums(int gid);
  void resetTransferSums(int gid);
};

