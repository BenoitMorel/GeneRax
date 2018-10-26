#ifndef JOINTSEARCH_UNDATEDDLMODEL_HPP_
#define JOINTSEARCH_UNDATEDDLMODEL_HPP_

#include <ale/containers/GeneMap.h>
#include <likelihoods/LibpllEvaluation.hpp>

using namespace std;

class UndatedDLModel {
public:
  int speciesNodesCount;

  // model
  vector<double> PD; // Duplication probability, per branch
  vector<double> PL; // Loss probability, per branch
  vector<double> PS; // Speciation probability, per branch
  const double O_R; // what is this?

  // SPECIES
  shared_ptr<bpp::PhyloTree> S;//Species tree
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
  void setRates(double dupRate, double lossRate);
  void setSpeciesTree(pll_rtree_t *geneTree);
  void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  double pun(shared_ptr<pllmod_treeinfo_t> treeinfo);
  

  void setMap(const GeneMap<string, string> &map);

  UndatedDLModel();
  ~UndatedDLModel();


  void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds);
  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  pll_unode_t *computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  void getRoots(pllmod_treeinfo_t &treeinfo, vector<pll_unode_t *> &roots);

};

#endif

