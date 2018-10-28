//all code by Szollosi GJ et al.; ssolo@elte.hu; GNU GPL 3.0;
#pragma once

#define ALE_VERSION "0.4"

#include <iostream>
#include <Bpp/Phyl/BipartitionTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTools.h>
#include <set>
#include <map>
#include<unordered_map>
#include <boost/progress.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>


#include "pairHasher.h"

//#include <boost/mpi.hpp>

typedef bpp::TreeTemplate<bpp::Node> tree_type;
//typedef long double scalar_type;
typedef  long double scalar_type;
typedef std::pair<bpp::Node*, bpp::Node*> dedge_type;
//typedef std::pair<long int,std::set < long int > > dip_type;


/****************************************************************************
 // approx_posterior class.
 // This class contains the description of a posterior on phylogenetic trees
 // by its clade probabilities and conditional clade probabilities.
 //Lexicon:
 // leaf set = bipartition = clade
 *****************************************************************************/

class approx_posterior
{ 
  
 public:
  //must load
  scalar_type observations;                    //Number of trees observed to build the approx_posterior object.


  //no need to load
  std::string constructor_string;              //string representing the tree in Newick format
  scalar_type alpha;
  scalar_type beta;
  std::map <std::string,int> tree_counts;

  ~approx_posterior()
    {
      leaf_ids.clear();
      id_leaves.clear();
      set_ids.clear(); 
      id_sets.clear();
      Bip_counts.clear();
      for (std::vector < std::unordered_map< std::pair<long int, long int> ,scalar_type> >::iterator it=Dip_counts.begin();it!=Dip_counts.end();it++)    
	(*it).clear();
      Dip_counts.clear();
      Gamma_s.clear();
      tree_bipstrings.clear();
      bipstring_trees.clear();
      set_sizes.clear();
      for (std::map <int, std::vector <long int > > :: iterator it = size_ordered_bips.begin(); it != size_ordered_bips.end(); it++)
	(*it).second.clear();
      size_ordered_bips.clear();
      Bip_bls.clear();
    }

  approx_posterior();                                                               //Does nothing. Formal constructor must be followed by load_state.
  approx_posterior(std::string tree);                                               //Constructs a basic instance by calling construct.
  void construct(std::string tree_string);                                          //Constructs a basic instance.

  void save_state(std::string fname) ;                                               //Writes the object to a file.
  void load_state(std::string fname);

  void observation(std::vector<std::string> trees,bool count_topologies=false, scalar_type weight=1.0);     //Given a vector of trees, fills an approx_posterior object by recursively calling decompose, and then doing some more counts.
  scalar_type p(std::string tree) const;                                                  //Computes the probability of a string tree. Calls recompose on the tree, and then uses the map returned by recompose to compute the probability of the whole tree.
  scalar_type nbipp(std::string tree) const;                                              //Computes the proportion of bipartitions already in the approx_posterior object that are present in the tree
  scalar_type binomial(int n,int m) const;                                                //Computes the binomial coefficient.
  scalar_type trinomial(int n1,int n2,int n3) const;                                      //Computes the multinomial coefficient for 3 elements.

  std::pair<std::string,scalar_type> mpp_tree() const;                                    //Returns the maximum a posteriori tree that can be amalgamated from the approx_prior object. Uses a double-recursive traversal of all bipartitions.
  std::string mpp_backtrack(long int g_id, std::map<long int, scalar_type > * qmpp) const;//Recursive function that, given a bipartition id and a map associating bipartition ids to their maximum a posteriori value, builds the maximum a posteriori tree, complete with (average) branch lengths.
  std::string random_tree() const;                                                        //Function that returns a random tree with unit branch lengths. Calls random_split.
	std::vector<std::string>  all_trees() const {return all_trees(Gamma);};                //del-loc. Builds all rooted trees that can be built based on the complete set of leaves.


  //no need to load
	std::set<int> Gamma_s;                                                             //del-loc. Set containing all leaf ids. ~Clade of all leaves in the tree.
    
    boost::dynamic_bitset<> Gamma;                                                                      //bit vector with all bits to 1 (all species present)
	int Gamma_size;                                                                  //Number of leaves.
   // size_t nbint;                                                                    //Number of ints in a bitvector used to store a partition
   // int filterForXor;                                                                //Integer used for xor operations to leave the unused bits of the last integer of the bitvector to 0
	scalar_type K_Gamma;                                                             //number of bipartitions of Gamma.
  long double N_Gamma;                                                               //number of unrooted trees on Gamma_size leaves
  std::string name_separator;                                                        //Character used between leaf names when writing the content of a leaf set.
  std::map <std::string,std::string> tree_bipstrings;                                //del-loc. Map between tree string and string containing all bipartitions in the tree.
  std::map <std::string,std::string> bipstring_trees;                                //del-loc. Dual from above. Map between string containing all bipartitions in the tree and tree string.
  std::map <long int,int> set_sizes;                                                 //del-loc. Map between a bipartition id and the sizes of the corresponding leaf set.
  std::map <int, std::vector <long int > > size_ordered_bips;                        //del-loc. Map between bipartition size, and the ids of all bipartitions of this size.

  //must load
  long int last_leafset_id;                                                          //Total number of sets of leaves (=bipartitions) observed in the posterior.
  std::map<std::string,int> leaf_ids;                                                //del-loc. Map between species name and leaf id. Leaf ids go from 1 to Gamma_size.
  std::map<int,std::string> id_leaves;                                               //del-loc. Map between leaf id and species name. Dual from above.
  std::map <long int,scalar_type> Bip_counts;                                        //del-loc. For each bipartition, gives the number of times it was observed.
  std::map <long int,scalar_type> Bip_bls;                                           //del-loc. Sum of the branch lengths associated to the bipartitions.

  //VECTORIZED BELOW  std::map <long int, std::map< std::set<long int>,scalar_type> > Dip_counts;        //del-loc. Contains the frequency of triplets: mother clade and its two daughter clades. Map between the bipartition id of the mother clade and another map containing a set of bipartition ids (couldn't it be just a pair, in the case of bifurcating trees?) and the frequency of the associated triplet of bipartitions. 

  std::vector < std::unordered_map< std::pair<long int, long int>,scalar_type> > Dip_counts;  
 // std::map <std::set <int>,long int>  set_ids;                                       //del-loc. Map between a set of leaf ids and the corresponding bipartition index.
//    std::map< long int, std::set <int> > id_sets;                                      //del-loc. Dual from above. Map between a bipartition index and the corresponding leaf ids.

    std::map < boost::dynamic_bitset<>,long int>  set_ids;                                   //del-loc. Map between a bit vector of leaf ids and the corresponding bipartition index.
    std::map< long int, boost::dynamic_bitset<> > id_sets;                                      //del-loc. Dual from above. Map between a bipartition index and the corresponding leaf id bit vector.

    
  //nuisance vars
  boost::timer * t;

  //algorithmic
  void decompose(std::string G_string,std::set<int> * bip_ids=NULL , scalar_type weight=1.0);                //Parses a tree in string format and updates the approx_prior object accordingly (notably updates the Bip_bls, Bip_counts, Dip_counts, and set_ids + id_sets through set2id)
  std::map < boost::dynamic_bitset<> ,scalar_type> recompose(std::string G_string) const;              //For a given input tree string, returns a map between all sets of leaves contained in the tree and their corresponding conditional clade probability.
  void register_leafset(std::string);
  long int set2id( boost::dynamic_bitset<>  leaf_set) ;                                           //If the set exists, returns the set id, otherwise creates a new set id for this set and returns it.

  //numeric
  scalar_type Bi(int n2) const;                                                            //Returns the total number of binary tree topologies possible given a fixed bipartition between n2 leaves on one side and Gamma_size-n2 leaves on the other side.
  scalar_type Tri(int n2,int n3) const;                                                    //Returns the total number of binary tree topologies possible given a fixed trifurcation between n2 leaves in one clade, n3 in another, and Gamma_size-n2-n3 leaves in the last one.

  scalar_type p_dip(long int g_id,long int gp_id,long int gpp_id) const;                   //Probability of a trifurcation given by the ids of the clades.
  //scalar_type p_dip(std::set<int> gamma,std::set<int> gammap,std::set<int> gammapp); //Probability of a trifurcation given by the leaf sets of the clades.
    scalar_type p_dip(boost::dynamic_bitset<> gamma, boost::dynamic_bitset<> gammap, boost::dynamic_bitset<> gammapp) const; //Probability of a trifurcation given by the leaf sets of the clades.
    
  scalar_type p_bip(long int g_id) const;                                                  //Probability of a bipartition given by its id. Uses the correction term alpha.
//  scalar_type p_bip(std::set<int> gamma);                                            //Probability of a bipartition given by its leaf set.
    scalar_type p_bip(boost::dynamic_bitset<> gamma) const;                                            //Probability of a bipartition given by its leaf set.


  //nuisance
  std::string set2name( boost::dynamic_bitset<> leaf_set) const;                                      //Prints the leaf names of leaves contained in a leaf set
  std::string random_split( boost::dynamic_bitset<> gamma) const;                                     //Recursive function that returns a random subtree given a leaf set as input and given the approx_posterior object. Can return clades never observed in the posterior sample.
  std::vector<std::string>  all_trees( boost::dynamic_bitset<> gamma) const;                        //del-loc. Builds all rooted trees that can be built with leaf set gamma.
  scalar_type count_trees() const;                                                         //Counts trees that can be amalgamated with the approx_posterior object with the complete leaf set, without actually building these trees.
  scalar_type count_trees(long int g_id) const;                                            //Counts trees that can be amalgamated with the leaf set with id g_id, without actually building these trees.

  scalar_type count_all_trees( boost::dynamic_bitset<> gamma) const;                                 //Counts all trees that can be built with the leaf set gamma, without actually building these trees.
	void setAlpha ( scalar_type a ) ;                                                 //Set the value for the alpha parameter used for normalizing counts
	void setBeta ( scalar_type b );                                                  //Set the value for the beta parameter used for normalizing counts
	std::vector < std::string > getLeafNames() const;										 //get a vector containing all leaf names in the ale.
	void computeOrderedVectorOfClades (std::vector <long int>&  ids, std::vector <long int>& id_sizes); //fills the ids and id_sizes maps, ids of clades ordered by their size.

private:
  ;
};



//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
//"Given a set S, the power set (or powerset) of S, written P(S), or 2S, is the set of all subsets of S."
template<typename Set> std::set<Set> powerset(const Set& s, size_t n);
//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
template<typename Set> std::set<Set> powerset(const Set& s);
