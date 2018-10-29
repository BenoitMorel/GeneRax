#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>

using namespace std;

class UndatedDTLModel {
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
  
  // SPECIES libpll
  vector <pll_rnode_t *> speciesNodes;
  
  // CLVs
  vector<vector<double> > uq;
  vector<double> ll; 
  
  map<string, string> geneNameToSpeciesName;
  map<string, int> speciesNameToId;
 

  vector<int> geneIds;
  vector<int> geneToSpecies;
  
  pll_unode_t *geneRoot;
public:
  void setRates(double dupRate, double lossRate, double transferRate);
  void setRoot(pll_unode_t * root) {geneRoot = root;}
  pll_unode_t *getRoot() {return geneRoot;}
  void setSpeciesTree(pll_rtree_t *geneTree);
  void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  double pun(shared_ptr<pllmod_treeinfo_t> treeinfo);
  

  void setGeneSpeciesMap(const GeneSpeciesMapping &map);

  UndatedDTLModel();
  ~UndatedDTLModel();


  void getIdsPostOrderRec(pll_unode_t *node, 
    vector<bool> &marked,
    vector<int> &nodeIds);
  void fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) ;
  void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds);
  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  pll_unode_t *computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  void getRoots(pllmod_treeinfo_t &treeinfo, vector<pll_unode_t *> &roots);

};


