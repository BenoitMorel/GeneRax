#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>
#include <vector>

using namespace std; 

class DatedDLModel: public AbstractReconciliationModel {
public:
  virtual ~DatedDLModel() {};
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  virtual void setSpeciesTree(pll_rtree_t *specieseTree);
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual void setGeneSpeciesMap(const GeneSpeciesMapping &map);
  virtual void setRoot(pll_unode_t * root);
  virtual pll_unode_t *getRoot();
  virtual double computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);

private:

  struct DDL_CLV {
    // probability of this directed gene
    // to be rooted at speciesId on subdivision 
    // subdivisionId
    // clv[speciesId][subdivisionId]
    vector<vector<double> > clv;
  };

  pll_unode_t *geneRoot;

  // branch subdivision lengths
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


private:
  void buildSubdivisionsRec(pll_rnode_t *node);

};


