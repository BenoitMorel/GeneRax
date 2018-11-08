#pragma once

#include <likelihoods/ale/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>

using namespace std;

class UndatedDTLModel: public AbstractReconciliationModel {
public:
  int speciesNodesCount;

  // model
  vector<double> PD; // Duplication probability, per branch
  vector<double> PL; // Loss probability, per branch
  vector<double> PT; // Transfer probability, per branch
  vector<double> PS; // Speciation probability, per branch
  const double O_R; // what is this?

  // SPECIES
  vector<int> daughter;
  vector<int> son;
  vector<double> uE; // Probability for a gene to become extinct on each brance
  vector<double> ll; 
  
  // SPECIES libpll
  vector <pll_rnode_t *> speciesNodes;
  
  // CLVs
  vector<vector<double> > uq;
  // for each gene node, average (over the species nodes)
  // transfer even probability
  vector<double> mPTuq;  
  
  // average (over the species nodes)
  // extinction probability
  double mPTE;

  map<string, string> geneNameToSpeciesName;
  map<string, int> speciesNameToId;
 

  vector<int> geneIds;
  vector<int> geneToSpecies;
  
  pll_unode_t *geneRoot;
public:
  UndatedDTLModel();
  virtual ~UndatedDTLModel();
  
  // unherited from parents
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  virtual void setSpeciesTree(pll_rtree_t *geneTree);
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual void setGeneSpeciesMap(const GeneSpeciesMapping &map);
  virtual void setRoot(pll_unode_t * root) {geneRoot = root;}
  virtual pll_unode_t *getRoot() {return geneRoot;}
  virtual double computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);

  void getIdsPostOrderRec(pll_unode_t *node, 
    vector<bool> &marked,
    vector<int> &nodeIds);
  void fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) ;
  void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds);
  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  void updateCLV(pll_unode_t *geneNode);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  pll_unode_t *computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  void getRoots(pllmod_treeinfo_t &treeinfo, vector<pll_unode_t *> &roots);

};


