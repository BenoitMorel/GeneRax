//all code by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
//modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "ALE.h"

#include <ale/containers/GeneMap.h>

using namespace std;

class exODT_DL_model {
public:
  int alpha;
  int last_branch;
  int last_rank;

  vector<int> daughter;
  vector<int> son;
  map<int, string> extant_species;                   //del-loc. Map between leaf id (0 to # of leaves) and leaf name.

  shared_ptr<tree_type> S;                                             //Species tree
  map<shared_ptr<bpp::PhyloNode>, int> node_ids;                         //Map between node and its id.


  map<string, shared_ptr<bpp::PhyloNode>> name_node;
  map<shared_ptr<bpp::PhyloNode>, string> node_name;

  vector<scalar_type> uE; // Probability for a gene to become extinct on each branch
  vector<vector<scalar_type> > uq;

  scalar_type PD; // Duplication probability, per branch
  scalar_type PL; // Loss probability, per branch
  scalar_type PS; // Speciation probability, per branch
  int speciesLastLeaf;

  map<long int, string> gid_sps;                                     //del-loc. Map between clade id (from the approx_posterior object) and species included in that clade.
  GeneMap<string, string> speciesGeneMap;                                     // Map between gene and species names.
  map<string, scalar_type> fraction_missing;
  vector<long int> g_ids;
  vector<long int> g_id_sizes;
  map<long int, int> g_id2i;

  scalar_type O_R;
  scalar_type delta;
  scalar_type lambda;


public:
  void setRates( 
    scalar_type OR,
    scalar_type dupRate,
    scalar_type lossRate) {
    O_R = OR;
    delta = dupRate;
    lambda = lossRate;
  }
  void construct_undated(const string &Sstring,
                         const string &fractionMissingFile = ""); //Constructs an object given a species tree and file containing fractions of missing genes per species.
  void calculate_undatedEs();

  void step_one(shared_ptr<approx_posterior> ale);
  void gene_species_mapping(shared_ptr<approx_posterior> ale);
  void inner_loop(bool g_is_a_leaf,
                  long int &g_id,
                  vector<int> &gp_is,
                  vector<long int> &gpp_is,
                  int i
                  );
  scalar_type pun(shared_ptr<approx_posterior> ale, bool verbose = false);


  //implemented in exODT.cpp
  void setMap(const GeneMap<string, string> &map_) { speciesGeneMap = map_; }

  exODT_DL_model();
  ~exODT_DL_model();
  void reset() {
    uE.clear();
    uq.clear();
    gid_sps.clear();
    fraction_missing.clear();
    g_ids.clear();
    g_id_sizes.clear();
    g_id2i.clear();
  }

};
