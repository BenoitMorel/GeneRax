//all code by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
//modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "ALE.h"

#include <ale/containers/GeneMap.h>

//using namespace std;
struct step {
  int e;
  int ep;
  int epp;
  scalar_type t;
  int rank;
  long int g_id;
  long int gp_id;
  long int gpp_id;
  std::string event;
};


/****************************************************************************
 // exODT_model class.
 // This class contains the description of a time-sliced species tree,
 // with parameters of duplication, transfer and loss, and a "population size",
 // and an approx_posterior object. This object can compute the probability of
 // an approx_posterior given the species tree and values for the ODTL parameters.
 // Contains lots of maps to map time slices to nodes, nodes to time slices,
 // nodes to ages...
 // Each time slice ends at a speciation node, with a particular age.
 // Each time slice has a rank, with the most recent one with rank 0,
 // and the older one with rank number_of_species.
 *****************************************************************************/
class exODT_model {
public:
  std::map<std::string, scalar_type> scalar_parameter;//del_loc
  std::map<std::string, std::vector<scalar_type> > vector_parameter;//del_loc
  std::map<std::string, std::string> string_parameter;//del_loc

  int alpha;
  int last_branch;
  int last_rank;
  std::shared_ptr<approx_posterior> ale_pointer;                            //Pointer to an approx_posterior object on which dynamic programming is performed in p for instance.

  bool transfers;

  std::vector<int> daughter;
  std::vector<int> son;

  std::map<int, std::string> extant_species;                   //del-loc. Map between leaf id (0 to # of leaves) and leaf name.


  std::shared_ptr<tree_type> S;                                             //Species tree
  std::map<std::shared_ptr<bpp::PhyloNode>, int> node_ids;                         //Map between node and its id.

  std::map<int, scalar_type> t_begin;                        //del-loc. Map between the id of a node, and the beginning of the branch that leads to it.
  std::map<int, scalar_type> t_end;                          //del-loc. Map between the id of a node, and the end of the time slice it defines, corresponding to the age of this node.

  std::map<int, std::vector<int> > time_slices;              //del-loc. Map between rank of time slice and indices of branches going through it. Terminating branch is last in the vector.
  std::map<int, std::vector<int> > branch_slices;            //del-loc. Map between a branch and all the time slices it traverses.
  std::map<int, std::vector<scalar_type> > time_slice_times; //del-loc. Map between rank of time slice and all the end times of the sub-slices inside this time slice.
  std::map<int, scalar_type> time_slice_begins;            //del-loc. Map between rank of time slice and begin time of this time slice.

  //Variables used for computing.
  std::map<int, std::map<scalar_type, scalar_type> > Ee;                       //del-loc. Probability (scalar value) that a gene present at a given time slice (whose rank is the int key) at time the first scalar key is getting extinct before reaching extant species.
  std::map<std::string, std::shared_ptr<bpp::PhyloNode>> name_node;
  std::map<std::shared_ptr<bpp::PhyloNode>, std::string> node_name;
  std::map<std::string, std::map<std::string, int> > ancestral_names;
  std::vector<std::vector<int> > ancestral;
  std::vector<std::vector<int> > ancestors; // contains the ancestors of a given branch; useful to forbid transfers to them.

  std::vector<scalar_type> fm; // Fraction of genes missing at the tips
  std::vector<scalar_type> uE; // Probability for a gene to become extinct on each branch
  scalar_type mPTE; // Mean probability across all branches for a gene to be transferred to branch h and then become extinct on that branch h
  std::vector<scalar_type> mPTE_ancestral_correction; // branch-wise adjustments of mPTE to obtain the branch-wise probability for a gene to be transferred to branch h and then become extinct on that branch h by doing mPTE - mPTE_ancestral_correction[e]. These branch-wise corrections are here to forbid transfers to ancestors of a branch.
  int root_i;
  std::vector<std::vector<scalar_type> > uq;
  std::vector<scalar_type> mPTuq;
  std::vector<std::vector<scalar_type> > mPTuq_ancestral_correction;

  std::vector<scalar_type> PD; // Duplication probability, per branch
  std::vector<scalar_type> PT; // Transfer probability, per branch
  std::vector<scalar_type> PL; // Loss probability, per branch
  std::vector<scalar_type> PS; // Speciation probability, per branch
  int last_leaf;
  std::map<int, std::map<scalar_type, scalar_type> > Ge;                       //del-loc. Probability (scalar value) that a gene present at a given time slice (whose rank is the int key) actually reaches extant species.
  std::map<long int, std::map<scalar_type, std::map<int, scalar_type> > > q; //del-loc. Map between clade id (from the approx_posterior object) and a map between the time of a subslice and a map between branch id and probability of the clade given the ODTL model.

  std::vector<std::vector<std::vector<std::map<int, scalar_type> > > > qvec;// NO del-loc !!

  std::map<long int, std::map<scalar_type, std::map<int, step> > > q_step;   //del-loc
  std::map<long int, std::string> gid_sps;                                     //del-loc. Map between clade id (from the approx_posterior object) and species included in that clade.
  GeneMap<std::string, std::string> speciesGeneMap;                                     // Map between gene and species names.

  std::map<std::string, scalar_type> MLRec_events;                            //del-loc
  std::map<std::string, std::vector<scalar_type> > branch_counts;             //del-loc
  std::vector<std::string> Ttokens;                                           //del-loc

  std::map<long int, std::vector<std::string> > gid_events;             //del-loc
  std::map<long int, std::vector<scalar_type> > gid_times;             //del-loc
  std::map<long int, std::vector<int> > gid_branches;             //del-loc
  std::map<long int, std::vector<long int> > gid_gidp;             //del-loc
  std::map<long int, std::vector<long int> > gid_gidpp;             //del-loc
  std::map<std::string, scalar_type> fraction_missing;

  void construct_undated(const std::string &Sstring,
                         const std::string &fractionMissingFile = ""); //Constructs an object given a species tree and file containing fractions of missing genes per species.
  void calculate_undatedEs();

  void clear_all();
  void step_one(std::shared_ptr<approx_posterior> ale);
  void gene_secies_mapping(std::shared_ptr<approx_posterior> ale);
  void inner_loop(std::shared_ptr<approx_posterior> ale, bool is_a_leaf,
                  long int &g_id,
                  std::vector<int> &gp_is,
                  std::vector<long int> &gpp_is,
                  std::vector<scalar_type> &p_part,
                  int i,
                  scalar_type &new_mPTuq);
  scalar_type pun(std::shared_ptr<approx_posterior> ale, bool verbose = false);
  //std::string feSPR(int e, int f);
  //std::vector<std::string> NNIs(int e);

  std::vector<long int> g_ids;
  std::vector<long int> g_id_sizes;
  std::map<long int, int> g_id2i;

  //implemented in exODT.cpp
  void setMap(const GeneMap<std::string, std::string> &map_) { speciesGeneMap = map_; }

  //void construct(const std::string& Sstring, const scalar_type& N=1e6, const std::string& fractionMissingFile=""); //Constructs an object given a species tree, population size and file containing fractions of missing genes per species.
  exODT_model();
  ~exODT_model();
  void set_model_parameter(std::string name, std::string value);    //Sets the value of a string parameter.
  void set_model_parameter(std::string, scalar_type);               //Sets the value of a scalar parameter.

};
