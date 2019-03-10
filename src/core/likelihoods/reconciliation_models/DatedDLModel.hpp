#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <vector>
#include <maths/ScaledValue.hpp>

using namespace std; 

class DatedDLModel: public AbstractReconciliationModel {
public:
  DatedDLModel();
  virtual ~DatedDLModel() {};
  // overloaded from parent
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  // overloaded from parent
  virtual bool implementsTransfers() {return false;}
protected:
  // overload from parent
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // overloaded from parent
  virtual void setSpeciesTree(pll_rtree_t *speciesTree);
  // overloaded from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overloaded from parent
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot);
  // overloaded from parent
  virtual ScaledValue getRootLikelihood(pll_unode_t *root) const;
  virtual ScaledValue getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {return getRecProba(root->node_index + _maxGeneId + 1, speciesRoot->node_index);}
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false) {};

private:
  double dupRate_;
  double lossRate_;
  double diffRates_;

  struct DDL_CLV {
    // probability of this directed gene
    // to be rooted at speciesId on subdivision 
    // subdivisionId
    // clv[speciesId][subdivisionId]
    vector<vector<ScaledValue> > clv;
  };


  // branch subdivision lengths
  // by convention, the first subdvision of each species 
  // has a 0 length, because it corresponds to the discontinuity
  // between this branch and its sons
  // then the subdivision go up to the parent node
  //
  // branchSubdivisions[speciesId][subdivisionId]
  vector<vector<double > > branchSubdivisions_;

  // gene extinction probabilities
  // extinctionProba_[speciesId][subdivisionId]
  vector<vector<double> > extinctionProba_;

  // single gene propagator probabilities
  // propagationProba_[speciesId][subdivisionId]
  vector<vector<double> > propagationProba_;

  // per gene intermediate result 
  // (conditional likelihood vectors, per analogy with libpll)
  // clvs_[geneId]
  vector<DDL_CLV> clvs_;
  DDL_CLV virtualCLV_;

  // probability that an extant gene is sampled
  // set to 1.0 fpr now
  double probaGeneSampled_;
  
  // map between extant genes and extant species
  vector<int> geneToSpecies;
private:
  
  void computeExtinctionProbas(pll_rtree_t *speciesTree);
  double propagateExtinctionProba(double initialProba, double branchLength); 
  void computePropagationProbas(pll_rtree_t *speciesTree);
  double propagatePropagationProba(double initialProba, double branchLength); 
  void computeCLVCell(pll_unode_t *geneNode, pll_rnode_t *speciesNode, vector<ScaledValue> &speciesCell, bool isVirtualRoot);
  void computeCLV(pll_unode_t *geneNode, pll_rnode_t *speciesNode, DDL_CLV *clv, bool isVirtualRoot = false);
  ScaledValue computeRecProbaInterBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode, bool isVirtualRoot);
  ScaledValue computeRecProbaIntraBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode, int subdivision, bool isVirtualRoot);
  const ScaledValue &getRecProba(int geneId, int speciesId) {
    return clvs_[geneId].clv[speciesId].back();
  }
  const ScaledValue &getRecProba(int geneId, int speciesId, int subdivision) {
    return clvs_[geneId].clv[speciesId][subdivision];
  }
  double getExtProba(int speciesId);
};

