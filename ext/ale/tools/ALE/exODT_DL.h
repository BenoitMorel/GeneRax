//all code by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
//modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "ALE.h"

#include <ale/containers/GeneMap.h>

class exODT_DL_model {
public:
  std::map<std::string, scalar_type> scalar_parameter;//del_loc
  std::map<std::string, std::vector<scalar_type> > vector_parameter;//del_loc
  std::map<std::string, std::string> string_parameter;//del_loc

  int alpha;
  int last_branch;
  int last_rank;

  std::vector<int> daughter;
  std::vector<int> son;
  std::map<int, std::string> extant_species;                   //del-loc. Map between leaf id (0 to # of leaves) and leaf name.

  std::shared_ptr<tree_type> S;                                             //Species tree
  std::map<std::shared_ptr<bpp::PhyloNode>, int> node_ids;                         //Map between node and its id.


  std::map<std::string, std::shared_ptr<bpp::PhyloNode>> name_node;
  std::map<std::shared_ptr<bpp::PhyloNode>, std::string> node_name;

  std::vector<scalar_type> uE; // Probability for a gene to become extinct on each branch
  std::vector<std::vector<scalar_type> > uq;

  scalar_type PD; // Duplication probability, per branch
  scalar_type PL; // Loss probability, per branch
  scalar_type PS; // Speciation probability, per branch
  int last_leaf;

  std::map<long int, std::string> gid_sps;                                     //del-loc. Map between clade id (from the approx_posterior object) and species included in that clade.
  GeneMap<std::string, std::string> speciesGeneMap;                                     // Map between gene and species names.
  std::map<std::string, scalar_type> fraction_missing;


public:
  void construct_undated(const std::string &Sstring,
                         const std::string &fractionMissingFile = ""); //Constructs an object given a species tree and file containing fractions of missing genes per species.
  void calculate_undatedEs();

  void step_one(std::shared_ptr<approx_posterior> ale);
  void gene_secies_mapping(std::shared_ptr<approx_posterior> ale);
  void inner_loop(std::shared_ptr<approx_posterior> ale, bool g_is_a_leaf,
                  long int &g_id,
                  std::vector<int> &gp_is,
                  std::vector<long int> &gpp_is,
                  int i
                  );
  scalar_type pun(std::shared_ptr<approx_posterior> ale, bool verbose = false);

  std::vector<long int> g_ids;
  std::vector<long int> g_id_sizes;
  std::map<long int, int> g_id2i;

  //implemented in exODT.cpp
  void setMap(const GeneMap<std::string, std::string> &map_) { speciesGeneMap = map_; }

  exODT_DL_model();
  ~exODT_DL_model();
  void set_model_parameter(std::string name, std::string value);    //Sets the value of a string parameter.
  void set_model_parameter(std::string, scalar_type);               //Sets the value of a scalar parameter.

};
