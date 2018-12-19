#ifndef JOINTSEARCH_UNDATEDDLMODEL_HPP_
#define JOINTSEARCH_UNDATEDDLMODEL_HPP_

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

  // model
  vector<double> PD; // Duplication probability, per branch
  vector<double> PL; // Loss probability, per branch
  vector<double> PS; // Speciation probability, per branch
  const double O_R; // what is this?

  // SPECIES
  vector<double> uE; // Probability for a gene to become extinct on each brance
  
  // CLVs
  vector<vector<ScaledValue> > uq;
  vector<ScaledValue> ll;
  
  vector<int> geneIds;

  bool allCLVInvalid;
  unordered_set<int> invalidCLVs;



  // repeats
  vector<unsigned long> repeatsId; // repeatsId[geneId]
  vector<vector<ScaledValue> > cache_;
public:
  UndatedDLModel();
  virtual ~UndatedDLModel();
  
  // unherited from parents
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  virtual double computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual void invalidateCLV(int nodeIndex);
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  

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

#endif

