//all code by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
//modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "ALE.h"

#include <ale/containers/GeneMap.h>

using namespace std;

class exODT_DL_model {
public:
  int last_branch;

  // model
  scalar_type PD; // Duplication probability, per branch
  scalar_type PL; // Loss probability, per branch
  scalar_type PS; // Speciation probability, per branch
  const scalar_type O_R; // what is this?

  // SPECIES
  shared_ptr<tree_type> S;                                             //Species tree
  vector<int> daughter;
  vector<int> son;
  map<int, string> extant_species;                   //del-loc. Map between leaf id (0 to # of leaves) and leaf name.
  int speciesLastLeaf;



  // CLVs
  vector<scalar_type> uE; // Probability for a gene to become extinct on each branch
  vector<vector<scalar_type> > uq;


  map<long int, string> gid_sps;                                     //del-loc. Map between clade id (from the approx_posterior object) and species included in that clade.
  GeneMap<string, string> speciesGeneMap;                                     // Map between gene and species names.
  vector<long int> g_ids;
  vector<long int> g_id_sizes;
  map<long int, int> g_id2i;



public:
  void setRates(scalar_type dupRate, scalar_type lossRate);
  void calculate_undatedEs();
  
  void construct_undated(const string &Sstring);
  

  void step_one(shared_ptr<approx_posterior> ale);
  void step_one_plop(shared_ptr<approx_posterior> ale);
  void gene_species_mapping(shared_ptr<approx_posterior> ale);
  void inner_loop(bool g_is_a_leaf,
                  long int &g_id,
                  vector<int> &gp_is,
                  vector<long int> &gpp_is,
                  int i
                  );
  scalar_type pun(shared_ptr<approx_posterior> ale);


  //implemented in exODT.cpp
  void setMap(const GeneMap<string, string> &map_) { speciesGeneMap = map_; }

  exODT_DL_model();
  ~exODT_DL_model();
  void reset() {
    uE.clear();
    uq.clear();
    gid_sps.clear();
    g_ids.clear();
    g_id_sizes.clear();
    g_id2i.clear();
  }

};
