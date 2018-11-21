#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>
#include <vector>

using namespace std; 

class DatedDLModel: public AbstractReconciliationModel {
public:
  DatedDLModel();
  virtual ~DatedDLModel() {};
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  virtual void setSpeciesTree(pll_rtree_t *speciesTree);
  virtual double computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);

private:

  double dupRate_;
  double lossRate_;
  double diffRates_;

  struct DDL_CLV {
    // probability of this directed gene
    // to be rooted at speciesId on subdivision 
    // subdivisionId
    // clv[speciesId][subdivisionId]
    vector<vector<double> > clv;
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
  void updateCLV(pll_unode_t *geneNode);
  double computeRecProbaInterBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode);
  double computeRecProbaIntraBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode, int subdivision);
  double getRecProba(int geneId, int speciesId);
  double getRecProba(int geneId, int speciesId, int subdivision);
  double getExtProba(int speciesId);
};


