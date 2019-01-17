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
* In addition, we forbid transfers to parent species
*/
class UndatedDTLModel: public AbstractReconciliationModel {
public:
  UndatedDTLModel();
  virtual ~UndatedDTLModel();
  
  // overloaded from parent
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);  
  virtual void invalidateCLV(int nodeIndex);
  virtual bool implementsTransfers() {return true;}
protected:
  // overloaded from parent
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual double computeLogLikelihoodInternal(shared_ptr<pllmod_treeinfo_t> treeinfo);
private:
  // model
  vector<double> _PD; // Duplication probability, per branch
  vector<double> _PL; // Loss probability, per branch
  vector<double> _PT; // Transfer probability, per branch
  vector<double> _PS; // Speciation probability, per branch

  // SPECIES
  vector<ScaledValue> _uE; // Probability for a gene to become extinct on each brance
  ScaledValue _transferExtinctionSum;
  vector<ScaledValue> _ancestralExctinctionCorrection;  

  // CLVs
  // _uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  vector<vector<ScaledValue> > _uq;
  vector<ScaledValue> _survivingTransferSums;
  vector<vector<ScaledValue> > _ancestralCorrection;

  // geme ids in postorder 
  vector<int> _geneIds;

  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  unordered_set<int> _invalidatedNodes;

  // is the CLV up to date?
  vector<bool> _isCLVUpdated;

  int _maxGeneId;

private:
  void getCLVsToUpdate(pllmod_treeinfo_t &treeinfo);
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot = false) const;
  void updateCLV(pll_unode_t *geneNode);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  void computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  void updateRoot(pllmod_treeinfo_t &treeinfo);
  double getSumLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  void invalidateAllCLVs();
  void updateCLVsRec(pll_unode_t *node);
  void markInvalidatedNodes(pllmod_treeinfo_t &treeinfo);
  void markInvalidatedNodesRec(pll_unode_t *node);
  void updateTransferSums(ScaledValue &transferExtinctionSum,
    vector<ScaledValue> &ancestralExtinctionCorrection,
    const vector<ScaledValue> &probabilities);
  void resetTransferSums(ScaledValue &transferSum,
    vector<ScaledValue> &ancestralCorrection);
  ScaledValue getCorrectedTransferExtinctionSum(int speciesNode) const;
  ScaledValue getCorrectedTransferSum(int geneId, int speciesId) const;
};

