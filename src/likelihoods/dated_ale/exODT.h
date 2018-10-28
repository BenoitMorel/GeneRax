//all code by Szollosi GJ et al.; ssolo@elte.hu; GNU GPL 3.0;
#include "ALE.h"
//using namespace std;
struct step
{
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
class exODT_model
{
 public:
  std::map <int,int> rank2label;

  std::map <std::string,scalar_type> scalar_parameter;//del_loc
  std::map <std::string,std::vector <scalar_type> > vector_parameter;//del_loc
  std::map <std::string,std::string> string_parameter;//del_loc

  int signal;
  std::string signal_string;

  int alpha;
  int last_branch;
  int last_rank;
  approx_posterior * ale_pointer;                            //Pointer to an approx_posterior object on which dynamic programming is performed in p for instance.

  //Runge-Krutta variables
  std::vector<scalar_type> Ee_y ;//del-loc

  scalar_type Ee_y_1;
  //  map<int,scalar_type> Ge_y;//del-loc
  std::vector<scalar_type> Ge_y ;//del-loc
  scalar_type Ge_y_1;

  std::vector<scalar_type> E_k1 ;
  std::vector<scalar_type> E_k2 ;
  std::vector<scalar_type> E_k3 ;
  std::vector<scalar_type> E_k4 ;//del-loc. Maps used for Runge-Kutta computations (4 stages).
  scalar_type E_k1_1, E_k2_1, E_k3_1, E_k4_1;

  std::vector<scalar_type> G_k1;
  std::vector<scalar_type> G_k2;
  std::vector<scalar_type> G_k3;
  std::vector<scalar_type> G_k4;//del-loc. Maps used for Runge-Kutta computations (4 stages).
  scalar_type G_k1_1, G_k2_1, G_k3_1, G_k4_1;


  std::map<int,int> father;                                  //del-loc. Map between node id and id of its father.
  std::map<int,std::vector<int> > daughters;                 //del-loc. Map between node id and ids of its daughters (-1 if node is a leaf).
  std::map<int,int> daughter;
  std::map<int,int> son;

  std::map<int,std::string> extant_species;                  //del-loc. Map between leaf id (0 to # of leaves) and leaf name.
  std::map <int,std::string> extant_taxa;                    // extant_taxa map for id-ing branches across trees


  std::map<int, scalar_type> branch_ts;                      //del-loc. Map between branch identified by the id of the node it ends at, and time of the time slice.
  std::map<int,int>rank_ids;                                 //del-loc. Map between rank of a time slice and the id of the node that terminates it.
  std::map<int,int>id_ranks;                                 //del-loc. Map between node id and rank of the time slice it terminates. Time slice at leaves has rank 0.

  tree_type * S;                                             //Species tree
  bpp::Node * S_root;                                        //Root of the species tree
  std::map<bpp::Node *,int>node_ids;                         //Map between node and its id.
  std::map<int,bpp::Node *>id_nodes;                         //Dual from the above, map between node id and node.

  std::map<int,scalar_type > t_begin;                        //del-loc. Map between the id of a node, and the beginning of the branch that leads to it.
  std::map<int,scalar_type > t_end;                          //del-loc. Map between the id of a node, and the end of the time slice it defines, corresponding to the age of this node.

  std::map<int,std::vector<int > > time_slices;              //del-loc. Map between rank of time slice and indices of branches going through it. Terminating branch is last in the vector.
  std::map<int,std::vector<int > > branch_slices;            //del-loc. Map between a branch and all the time slices it traverses.
  std::map<int,std::vector<scalar_type > > time_slice_times; //del-loc. Map between rank of time slice and all the end times of the sub-slices inside this time slice.
	std::map<int,scalar_type > time_slice_begins;            //del-loc. Map between rank of time slice and begin time of this time slice.

  //Variables used for computing.
  std::map<int,std::map <scalar_type,scalar_type> > Ee;                       //del-loc. Probability (scalar value) that a gene present at a given time slice (whose rank is the int key) at time the first scalar key is getting extinct before reaching extant species.
  std::map <std::string,bpp::Node *> name_node;
  std::map <bpp::Node *,std::string> node_name;
  std::map <std::string,std::map<std::string,int> > ancestral_names;
  std::map <int,std::map<int,int> > ancestral;
  std::vector < std::vector <int> > ancestors; // contains the ancestors of a given branch; useful to forbid transfers to them.

  std::vector<scalar_type> fm; // Fraction of genes missing at the tips
  std::vector<scalar_type> uE; // Probability for a gene to become extinct on each branch
  scalar_type mPTE; // Mean probability across all branches for a gene to be transferred to branch h and then become extinct on that branch h
  std::vector<scalar_type> mPTE_ancestral_correction; // branch-wise adjustments of mPTE to obtain the branch-wise probability for a gene to be transferred to branch h and then become extinct on that branch h by doing mPTE - mPTE_ancestral_correction[e]. These branch-wise corrections are here to forbid transfers to ancestors of a branch.
  int root_i;
  std::vector < std::vector <scalar_type> > uq;
  std::vector < scalar_type > mPTuq;
  std::vector < std::vector <scalar_type> > mPTuq_ancestral_correction;

  std::vector<scalar_type> PD; // Duplication probability, per branch
  std::vector<scalar_type> PT; // Transfer probability, per branch
  std::vector<scalar_type> PL; // Loss probability, per branch
  std::vector<scalar_type> PS; // Speciation probability, per branch
  int last_leaf;
  std::map<int,std::map <scalar_type,scalar_type> > Ge;                       //del-loc. Probability (scalar value) that a gene present at a given time slice (whose rank is the int key) actually reaches extant species.
  std::map<long int, std::map< scalar_type, std::map<int, scalar_type> > > q; //del-loc. Map between clade id (from the approx_posterior object) and a map between the time of a subslice and a map between branch id and probability of the clade given the ODTL model.

  std::vector< std::vector < std::vector < std::map<int, scalar_type> > > > qvec;// NO del-loc !!

  std::map<long int, std::map< scalar_type, std::map<int, step> > > q_step;   //del-loc
  std::map <long int,std::string> gid_sps;                                    //del-loc. Map between clade id (from the approx_posterior object) and species included in that clade.

  std::map <std::string,scalar_type> MLRec_events;                            //del-loc
  std::map<std::string, std::vector<scalar_type> > branch_counts;             //del-loc
  std::vector<std::string> Ttokens;                                           //del-loc

  std::map<long int, std::vector<std::string> > gid_events;             //del-loc
  std::map<long int, std::vector<scalar_type> > gid_times;             //del-loc
  std::map<long int, std::vector<int> > gid_branches;             //del-loc
  std::map<long int, std::vector<long int> > gid_gidp;             //del-loc
  std::map<long int, std::vector<long int> > gid_gidpp;             //del-loc
  std::map<std::string, scalar_type> fraction_missing;

  void construct_undated(const std::string& Sstring, const std::string& fractionMissingFile=""); //Constructs an object given a species tree and file containing fractions of missing genes per species.
  void calculate_undatedEs();
  scalar_type pun(approx_posterior *ale, bool verbose=false);
  std::string feSPR(int e, int f);
  std::vector<std::string> NNIs(int e);

  std::string sample_undated();
  std::string sample_undated(int e,int i,std::string last_event,std::string branch_string="");
  std::vector <long int>  g_ids;
  std::vector <long int>  g_id_sizes;
  std::map <long int,int> g_id2i;

  //implemented in exODT.cpp

  void construct(const std::string& Sstring, const scalar_type& N=1e6, const std::string& fractionMissingFile=""); //Constructs an object given a species tree, population size and file containing fractions of missing genes per species.
  exODT_model();
  ~exODT_model()
    {
      rank_ids.clear();
      id_ranks.clear();
      father.clear();
      for (std::map<int,std::vector<int> >::iterator it=daughters.begin();it!=daughters.end();it++)
	(*it).second.clear();
      daughters.clear();
      extant_species.clear();
      branch_ts.clear();
      rank_ids.clear();
      id_ranks.clear();
      t_begin.clear();
      t_end.clear();
      for (std::map<int,std::vector<int> >::iterator it=time_slices.begin();it!=time_slices.end();it++)
	(*it).second.clear();
      time_slices.clear();
      for (std::map<int,std::vector<int> >::iterator it=branch_slices.begin();it!=branch_slices.end();it++)
	(*it).second.clear();
      for (std::map<int,std::vector<scalar_type > >::iterator it=time_slice_times.begin();it!=time_slice_times.end();it++)
	(*it).second.clear();
      time_slice_times.clear();
      time_slice_begins.clear();
      scalar_parameter.clear();
      for (std::map <std::string,std::vector <scalar_type> >::iterator it=vector_parameter.begin();it!=vector_parameter.end();it++)//del_loc
	(*it).second.clear();
      vector_parameter.clear();
      string_parameter.clear();
      node_ids.clear();
      id_nodes.clear();
      delete S;
      for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ee.begin();it!=Ee.end();it++)//del_loc
	(*it).second.clear();
      Ee.clear();
      for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ge.begin();it!=Ge.end();it++)//del_loc
	(*it).second.clear();
      Ge.clear();
      Ee.clear();
      for (std::map<long int, std::map< scalar_type, std::map<int, scalar_type> > >::iterator it=q.begin();it!=q.end();it++)
	{
	  for ( std::map< scalar_type, std::map<int, scalar_type> >::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	    (*jt).second.clear();
	  (*it).second.clear();
	}
      q.clear();
      for (std::map<long int, std::map< scalar_type, std::map<int, step> > >::iterator it=q_step.begin();it!=q_step.end();it++)
	{
	  for ( std::map< scalar_type, std::map<int, step> >::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	    (*jt).second.clear();
	  (*it).second.clear();
	}
      q_step.clear();
      gid_sps.clear();
      MLRec_events.clear();
      for (std::map<std::string, std::vector<scalar_type> >::iterator it=branch_counts.begin();it!=branch_counts.end();it++)//del_loc
	(*it).second.clear();
      branch_counts.clear();

      for (std::map<long int, std::vector<std::string> >::iterator it=gid_events.begin();it!=gid_events.end();it++)//del_loc
	(*it).second.clear();
      gid_events.clear();

      for (std::map<long int, std::vector<scalar_type> >::iterator it=gid_times.begin();it!=gid_times.end();it++)//del_loc
	(*it).second.clear();
      gid_times.clear();

      for (std::map<long int, std::vector<int> >::iterator it=gid_branches.begin();it!=gid_branches.end();it++)//del_loc
	(*it).second.clear();
      gid_branches.clear();


      Ttokens.clear();
    };

  void set_model_parameter(std::string name,std::string value);    //Sets the value of a string parameter.
  void set_model_parameter(std::string,scalar_type);               //Sets the value of a scalar parameter.
  void set_model_parameter(std::string,std::vector<scalar_type>);  //Sets the value of a vector of scalars parameter.

  //implemented in model.cpp
  scalar_type p(approx_posterior *ale);                            //Computes the probability of an approx_posterior according to the species tree and parameter values.
  void calculate_EG();
  void calculate_EGb(); 										   //Fills Ee. Calculates extinction probabilities per branch and per time slice.

  //implemented in traceback.cpp
  std::pair<std::string,scalar_type> p_MLRec(approx_posterior *ale,bool lowmem=true);
  std::pair<std::string,scalar_type> traceback();
  std::string traceback(long int g_id,scalar_type t,scalar_type rank,int e,scalar_type branch_length,std::string branch_events,std::string transfer_token="");
  void register_O(int e);
  void register_D(int e);
  void register_Tto(int e);
  void register_Tfrom(int e);
  void register_L(int e);
  void register_S(int e);
  void register_Su(int e,std::string last_event);
  void register_T_to_from(int e,int f);
  void reset_T_to_from( );

  std::vector < std::vector<scalar_type> > T_to_from;

  void register_leaf(int e);
  void register_leafu(int e,std::string last_event);

  void register_Ttoken(std::string token);
  //implemented in traceback_lowmem.cpp - not done
  std::pair<std::string,scalar_type> p_MLRec_lowmem(approx_posterior *ale);
  std::string traceback_lowmem(long int g_id,scalar_type t,scalar_type rank,int e,scalar_type branch_length,std::string branch_events,std::string transfer_token="");
  //implemented in sample.cpp
  std::string sample(bool max_rec=false);
  std::string sample(bool S_node,long int g_id,int t_i,scalar_type rank,int e,scalar_type branch_length,std::string branch_events, std::string transfer_toke="",bool max_rec=false);


  void show_counts(std::string name, bool as_branch_length=true, bool per_copy=false);
  std::string counts_string(scalar_type samples=1);
  std::string counts_string_undated(scalar_type samples=1);

  void show_rates(std::string name);
  std::string gid_string(long int g_id);
  std::string vertical_string(long int g_id,std::string ancestral_string="",scalar_type t_0=-1);

 private:
  ;

};
