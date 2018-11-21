#ifndef JOINTSEARCH_UNDATEDDLMODEL_HPP_
#define JOINTSEARCH_UNDATEDDLMODEL_HPP_

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>

using namespace std;

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
  vector<vector<double> > uq;
  vector<double> ll; 
  
  vector<int> geneIds;
public:
  UndatedDLModel();
  virtual ~UndatedDLModel();
  
  // unherited from parents
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  virtual double computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  
private:
  void updateCLV(pll_unode_t *geneNode);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  pll_unode_t *computeLikelihoods(pllmod_treeinfo_t &treeinfo);

};

#endif

