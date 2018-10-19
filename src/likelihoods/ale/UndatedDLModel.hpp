#ifndef JOINTSEARCH_UNDATEDDLMODEL_HPP_
#define JOINTSEARCH_UNDATEDDLMODEL_HPP_

#include <ale/containers/GeneMap.h>
#include <likelihoods/LibpllEvaluation.hpp>

using namespace std;

class UndatedDLModel {
public:
  int last_branch;

  // model
  double PD; // Duplication probability, per branch
  double PL; // Loss probability, per branch
  double PS; // Speciation probability, per branch
  const double O_R; // what is this?

  // SPECIES
  shared_ptr<bpp::PhyloTree> S;//Species tree
  vector<int> daughter;
  vector<int> son;
  int speciesLastLeaf;
  vector<double> uE; // Probability for a gene to become extinct on each brance

  // CLVs
  vector<vector<double> > uq;
  vector<double> ll; 
  
  map<string, string> geneNameToSpeciesName;
  map<string, int> speciesNameToId;
 

  vector<int> geneIds;
  vector<int> inverseGeneIds;
  vector<int> geneToSpecies;
  
  /*
  vector<long int> g_ids;
  vector<long int> g_id_sizes;
  map<long int, int> g_id2i;
*/


public:
  void setRates(double dupRate, double lossRate);
  void construct_undated(const string &Sstring);
  double pun(shared_ptr<pllmod_treeinfo_t> treeinfo);
  

  void setMap(const GeneMap<string, string> &map);

  UndatedDLModel();
  ~UndatedDLModel();
  void reset() {
    uE.clear();
    uq.clear();
  }


  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  void computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  void getRoots(pllmod_treeinfo_t &treeinfo, vector<pll_unode_t *> &roots);

};

#endif

