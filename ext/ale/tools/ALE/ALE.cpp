// Modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "ALE.h"


#include <ale/tools/Utils.h>
#include <ale/tools/PhyloTreeToolBox.h>
#include <ale/tools/IO/IO.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/NumConstants.h>

using namespace std;
using namespace bpp;

#include <bitset>

//## aux. functions ##

//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
//"Given a set S, the power set (or powerset) of S, written P(S), or 2S, is the set of all subsets of S."
template<typename Set>
set<Set> powerset(const Set &s, size_t n) {
  typedef typename Set::const_iterator SetCIt;
  typedef typename set<Set>::const_iterator PowerSetCIt;
  set<Set> res;
  if (n > 0) {
    set<Set> ps = powerset(s, n - 1);
    for (PowerSetCIt ss = ps.begin(); ss != ps.end(); ss++)
      for (SetCIt el = s.begin(); el != s.end(); el++) {
        Set subset(*ss);
        subset.insert(*el);
        res.insert(subset);
      }
    res.insert(ps.begin(), ps.end());
  } else
    res.insert(Set());
  return res;
}

//from: http://rosettacode.org/wiki/Power_set#Recursive_version (accessed 12/13/11)
template<typename Set>
set<Set> powerset(const Set &s) {
  return powerset(s, s.size());
}



approx_posterior::approx_posterior(string tree_string) {
  constructor_string = tree_string;
  construct(constructor_string);
}


void approx_posterior::construct(string tree_string) {
  //  t=new boost::timer();

  last_leafset_id = 0;
  observations = 0;
  vector<string> leaves;
  // ? name_separator="+";

  // load leaf names.
  if (tree_string.substr(0, 1) != "(") {
    //boost::trim(tree_string);
    Utils::trim_str(tree_string);
    //boost::split(leaves, tree_string, boost::is_any_of(","), boost::token_compress_on);
    leaves = Utils::splitString(tree_string, ",");

  } else {
    std::shared_ptr<tree_type> tree(IO::newickToPhyloTree(tree_string, false));//del-loc

    leaves = tree->getAllLeavesNames();//del-loc
  }

  // give an id for each leaf.
  int id = 0;
  for (vector<string>::iterator it = leaves.begin(); it != leaves.end(); it++) {
    id++;
    string leaf_name = (*it);
    leaf_ids[leaf_name] = id;
    Gamma_s.insert(id);
    id_leaves[id] = leaf_name;
  }

  alpha = 0;
  beta = 0;
  Gamma_size = Gamma_s.size();
  Gamma = boost::dynamic_bitset<>(Gamma_size + 1);//new int[nbint];
  //All leaves are present in Gamma:
  for (auto i = 0; i < Gamma_size + 1; i++) {
    Gamma[i] = 1;
  }

  //maybe should use boost pow
  //number of bipartitions of Gamma
//  K_Gamma = pow(2., (int) Gamma_size - 1) - 1;
  //number of unrooted trees on Gamma_size leaves
  if (Gamma_size < 3)
    N_Gamma = 1;
  else
    N_Gamma = Utils::double_factorial<long double>(2 * Gamma_size - 5);

  //XX
  std::unordered_map<pair<long int, long int>, scalar_type> temp;
  while (leaves.size() + 1 > Dip_counts.size())
    Dip_counts.push_back(temp);
  //XX
  //del-locs
  leaves.clear();
  //delete tree;
}


long int approx_posterior::set2id(boost::dynamic_bitset<> leaf_set) {
  long int id = set_ids[leaf_set];
  if (!id) {
    last_leafset_id++;
    set_ids[leaf_set] = last_leafset_id;
    // TMP for debug
    //Dip_levels[leaf_set.size()].push_back(last_leafset_id);
    //id2name[last_leafset_id]=set2name(leaf_set);
    //VEC
    std::unordered_map<pair<long int, long int>, scalar_type> tmp;
    Dip_counts.push_back(tmp);
    //VEC
    id_sets[last_leafset_id] = leaf_set;
    Bip_bls[last_leafset_id] = 0;
    // std::cout << "notid: "<< last_leafset_id <<std::endl;
    return last_leafset_id;
  } else {
    // std::cout << "id: "<< id <<std::endl;
    return id;
  }
}


scalar_type approx_posterior::Bi(int n2) const {
  int n1 = Gamma_size - n2;
  if (n2 == 1 or n1 == 1)
    return Utils::double_factorial<scalar_type>(2 * Gamma_size - 5);
  n1 = max(2, n1);
  n2 = max(2, n2);
  return Utils::double_factorial<scalar_type>(2 * n1 - 3) *
         Utils::double_factorial<scalar_type>(2 * n2 - 3);
}


scalar_type approx_posterior::Tri(int n2, int n3) const {
  int n1 = Gamma_size - n2 - n3;
  n1 = max(2, n1);
  n2 = max(2, n2);
  n3 = max(2, n3);
  return Utils::double_factorial<scalar_type>(2 * n1 - 3) *
         Utils::double_factorial<scalar_type>(2 * n2 - 3) *
         Utils::double_factorial<scalar_type>(2 * n3 - 3);
}


scalar_type approx_posterior::p_bip(boost::dynamic_bitset<> gamma) const {
  if (Gamma_size < 4)
    return 1;
  long int g_id;
  if (set_ids.count(gamma))
    g_id = set_ids.at(gamma);
  else
    g_id = -10;
//    std::cout << "g_id: "<< g_id << " TO "<< p_bip(g_id) <<std::endl;
  return p_bip(g_id);
}


scalar_type approx_posterior::p_dip(boost::dynamic_bitset<> gamma, boost::dynamic_bitset<> gammap,
                                    boost::dynamic_bitset<> gammapp) const {
  if (Gamma_size < 4)
    return 1;
  long int g_id;
  if (set_ids.count(gamma))
    g_id = set_ids.at(gamma);
  else
    g_id = -10;

  long int gp_id;
  if (set_ids.count(gammap))
    gp_id = set_ids.at(gammap);
  else
    gp_id = -10;

  long int gpp_id;
  if (set_ids.count(gammapp))
    gpp_id = set_ids.at(gammapp);
  else
    gpp_id = -10;

  //std::cout << "g_id: "<< g_id << " AND "<< gp_id << " AND "<< gpp_id << " TO: "<< p_dip(g_id, gp_id, gpp_id) <<" " << p_dip( g_id, gpp_id, gp_id) <<std::endl;
  if (gpp_id > gp_id)
    return p_dip(g_id, gp_id, gpp_id);
  else
    return p_dip(g_id, gpp_id, gp_id);
}


scalar_type approx_posterior::p_bip(long int g_id) const {
  if (Gamma_size < 4)
    return 1;

  scalar_type Bip_count = 0;

  if (Bip_counts.count(g_id) == 0 or g_id == -10 or !g_id) {
    //never saw gamma in sample
    Bip_count = 0;
  } else {
    Bip_count = Bip_counts.at(g_id);
  }
  //if ( gamma.size()==1 or (int)gamma.size()==Gamma_size-1) Bip_count=observations;
  if (set_sizes.count(g_id) == 0 or g_id == -10) Bip_count = 0;
  else if (set_sizes.at(g_id) == 1 or set_sizes.at(g_id) == Gamma_size - 1) Bip_count = observations;

  if (alpha > 0)
    return Bip_count / (observations + alpha) + (alpha / N_Gamma * Bi(set_sizes.at(g_id))) / (observations + alpha);
  else
    return Bip_count / observations;
}


scalar_type approx_posterior::p_dip(long int g_id, long int gp_id, long int gpp_id) const {
  if (Gamma_size < 4)
    return 1;
  scalar_type beta_switch = 1;
  scalar_type Dip_count = 0, Bip_count = 0;
  if (Bip_counts.count(g_id) == 0 or g_id == -10 or !g_id) {
    //never saw gamma in sample
    beta_switch = 0.;
    Bip_count = 0;
    Dip_count = 0;
  } else {
    //set <long int> parts;
    //parts.insert(gp_id);
    //parts.insert(gpp_id);
    pair<long int, long int> parts;
    if (gpp_id > gp_id) {
      parts.first = gp_id;
      parts.second = gpp_id;
    } else {
      parts.first = gpp_id;
      parts.second = gp_id;
    }

    Bip_count = Bip_counts.at(g_id);
    //Dip_count=Dip_counts.at(g_id).at(parts);
    if (gp_id == -10 or gpp_id == -10 or Dip_counts.at(g_id).count(parts) == 0 or !gp_id or !gpp_id) {
      //never saw gammap-gammapp partition in sample
      Dip_count = 0;
    } else {
      Dip_count = Dip_counts.at(g_id).at(parts);
    }
  }
  if (set_sizes.count(g_id) == 0 or set_sizes.at(g_id) == 1 or set_sizes.at(g_id) == Gamma_size - 1)
    Bip_count = observations;
  // ? above is correct or not ?
  //cout << "Dip:"<<Dip_count << "  Bip:" << Bip_count << endl;
  if (alpha > 0 or beta > 0)
    return (Dip_count + (alpha / N_Gamma * Tri(set_sizes.at(gp_id), set_sizes.at(gpp_id))) +
            beta_switch * beta / (pow(2., set_sizes.at(g_id) - 1) - 1)) /
           (Bip_count + (alpha / N_Gamma * Bi(set_sizes.at(g_id))) + beta_switch * beta);
  else
    return Dip_count / Bip_count;
}


// an unrooted tree given by its Newick string (which can be rooted)
map<boost::dynamic_bitset<>, scalar_type> approx_posterior::recompose(string G_string) const {
  map<boost::dynamic_bitset<>, scalar_type> return_map;
  map<dedge_type, int> dedges;//del-loc
  map<Node, vector<Node> > neighbor;//del-loc

  std::shared_ptr<tree_type> G = IO::newickToPhyloTree(G_string, false);//del-loc

  //if (G->isRooted()) G->unroot();
  vector<Node> nodes = G->getAllNodes(); //del-loc

  //Find all directed edges
  for (vector<Node>::iterator it = nodes.begin(); it != nodes.end(); it++) {
    Node from = *it;
    if (G->hasFather(from)){//from->hasFather()) {
      Node father = G->getFather(from);//from->getFather();
      neighbor[from].push_back(father);
    }
    if (!G->isLeaf(from)) {
      vector<Node > sons = G->getSons(from);//from->getSons(); //del-loc
      for (vector<Node>::iterator it_sons = sons.begin(); it_sons != sons.end(); it_sons++)
        neighbor[from].push_back((*it_sons));
      sons.clear();
    }
    for (vector<Node>::iterator it_tos = neighbor[from].begin(); it_tos != neighbor[from].end(); it_tos++) {
      dedge_type dedge;
      dedge.first = from;
      dedge.second = *it_tos;
      dedges[dedge] = 0;
    }
  }
  nodes.clear();

  map<dedge_type, boost::dynamic_bitset<> > flat_names; //del-loc
  map<dedge_type, scalar_type> q; //del-loc
  // Name all leaves
  for (map<dedge_type, int>::iterator it = dedges.begin(); it != dedges.end(); it++) {

    dedge_type dedge = (*it).first;
    Node from = dedge.first;
    Node to = dedge.second;
    //visit dedges from a leaf as these can be named
    if (G->isLeaf(from)) {
      if (flat_names.find(dedge) == flat_names.end()) {
        //int* temp = new int[nbint];
        boost::dynamic_bitset<> temp(Gamma_size + 1);
        /*   for (auto i=0; i< nbint; ++i) { //Resetting all bits
               temp[i] = 0;
               //BipartitionTools::bit0( temp, static_cast<int>(i) );
           }*/
        flat_names[dedge] = temp;
      }
      //BipartitionTools::bit1(flat_names.at(dedge), static_cast<int>( leaf_ids.at(from->getName() ) ) );
      flat_names.at(dedge)[static_cast<int>( leaf_ids.at(from->getName()))] = 1;
      // flat_names[dedge].insert(leaf_ids[from->getName()]);
      q[dedge] = 1;
      return_map[flat_names[dedge]] = q[dedge];
      //mark named
      dedges[dedge] = -1;
      //proceed to dedges in next level -  dedges from cherries can now be named and at least one cherry must exist
      for (vector<Node>::iterator it_tos = neighbor[to].begin(); it_tos != neighbor[to].end(); it_tos++)
        if ((*it_tos) != from) {
          dedge_type dedge_out;
          dedge_out.first = to;
          dedge_out.second = *it_tos;
          dedges[dedge_out] += 1;
        }
    }
  }

  bool edges_left = false;
  for (map<dedge_type, int>::iterator it = dedges.begin(); it != dedges.end(); it++) {
    dedge_type dedge = (*it).first;
    if (dedges[dedge] != -1)
      edges_left = true;  //Couldn't we add a "break" here?
  }
  while (edges_left) {
    for (map<dedge_type, int>::iterator it = dedges.begin(); it != dedges.end(); it++) {
      dedge_type dedge = (*it).first;
      //Process edges that can be named
      if (dedges[dedge] == 2) {
        Node from = dedge.first;
        Node to = dedge.second;
        vector<dedge_type> dedges_in; //del-loc
        for (vector<Node>::iterator it_tos = neighbor[from].begin(); it_tos != neighbor[from].end(); it_tos++)
          if (*it_tos != to) {
            dedge_type dedge_in;
            dedge_in.first = *it_tos;
            dedge_in.second = from;
            dedges_in.push_back(dedge_in);
          }

        boost::dynamic_bitset<> leaf_set_in_1(Gamma_size + 1);
        leaf_set_in_1 = flat_names[dedges_in[0]];
        boost::dynamic_bitset<> leaf_set_in_2(Gamma_size + 1);
        leaf_set_in_2 = flat_names[dedges_in[1]];
        if (flat_names.find(dedge) == flat_names.end()) {
          //int* temp = new int[nbint];
          boost::dynamic_bitset<> temp(Gamma_size + 1);
          flat_names[dedge] = temp;
        }
        //  BipartitionTools::bitOr(flat_names[dedge], leaf_set_in_1, leaf_set_in_2, nbint);
        flat_names[dedge] = leaf_set_in_1 | leaf_set_in_2;

        q[dedge] = q[dedges_in[0]] * q[dedges_in[1]] * p_dip(flat_names[dedge], leaf_set_in_1, leaf_set_in_2);
        return_map[flat_names[dedge]] = q[dedge];

        //mark named
        dedges[dedge] = -1;
        //proceed to dedges in next level -  new dedges can now be named
        for (vector<Node>::iterator it_tos = neighbor[to].begin(); it_tos != neighbor[to].end(); it_tos++)
          if ((*it_tos) != from) {
            dedge_type dedge_out;
            dedge_out.first = to;
            dedge_out.second = *it_tos;
            dedges[dedge_out] += 1;
          }
        dedges_in.clear();
      }
    }
    edges_left = false;
    for (map<dedge_type, int>::iterator it = dedges.begin(); it != dedges.end(); it++) {
      dedge_type dedge = (*it).first;
      if (dedges[dedge] != -1)
        edges_left = true;
    }
  }

  //del-locs
  dedges.clear();
  for (map<Node, vector<Node> >::iterator it = neighbor.begin(); it != neighbor.end(); it++)
  (*it).second.clear();
  neighbor.clear();
  for (auto it = flat_names.begin(); it != flat_names.end(); it++) {
    for (auto i = 0; i < Gamma_size + 1; ++i) { //Resetting all bits
      (*it).second[i] = 0;
    }
  }
  flat_names.clear();
  q.clear();
  return return_map;
}


// an unrooted tree given by its Newick string (which can be rooted)
void approx_posterior::decompose(string G_string, set<int> *bip_ids, scalar_type weight) {
  //VEC
  std::unordered_map<pair<long int, long int>, scalar_type> tmp;
  Dip_counts.push_back(tmp);
  //VEC

  //vector <dip_type > return_dips;
  map<dedge_type, int> dedges;//del-loc
  map<Node, vector<Node> > neighbor;//del-loc

  std::shared_ptr<tree_type> G = IO::newickToPhyloTree(G_string, false);//del-loc

  if (G->isRooted()) {

    // unroot the tree
    auto root = G->getRoot();
    auto root_sons = G->getSons(root);
    if(root_sons.size() == 2) {
      G->rootAt(root_sons.front());
      for (std::size_t i = 1; i < root_sons.size(); i++) {
        std::shared_ptr<bpp::PhyloBranch> new_branch(new bpp::PhyloBranch());
        if(PhyloTreeToolBox::hasDistanceBetweenTwoNodes(root_sons.front(), root_sons.at(i), *G)) {
          double distance = PhyloTreeToolBox::getDistanceBetweenTwoNodes(root_sons.front(), root_sons.at(i), *G);
          new_branch->setLength(distance);
        }
        G->link(root_sons.front(), root_sons.at(i), new_branch);
        G->unlink(root, root_sons.at(i));
      }
      G->unlink(root_sons.front(), root);
      G->deleteNode(root);
    }
    if(G->getNumberOfSons(G->getRoot()) == 1) G->rootAt(G->getSons(G->getRoot()).front());

  }

  std::vector<Node> nodes = G->getAllNodes();

  //Constitute all directed edges in dedges.
  for (auto it = nodes.begin(); it != nodes.end(); it++) {
    Node from = *it;
    if (G->hasFather(from)){//from->hasFather()) {
      Node father = G->getFather(from);//from->getFather();
      neighbor[from].push_back(father);
    }
    if (not G->isLeaf(from)) {
      vector<Node> sons = G->getSons(from); //del-loc
      for (auto it_sons = sons.begin(); it_sons != sons.end(); it_sons++)
        neighbor[from].push_back((*it_sons));
      sons.clear();
    }
    for (auto it_tos = neighbor[from].begin(); it_tos != neighbor[from].end(); it_tos++) {
      dedge_type dedge;
      dedge.first = from;
      dedge.second = *it_tos;
      dedges[dedge] = 0;
    }
  }
  nodes.clear();

  map<dedge_type, boost::dynamic_bitset<> > flat_names; //del-loc
  // Name all leaves
  for (auto it = dedges.begin(); it != dedges.end(); it++) {

    dedge_type dedge = (*it).first;
    Node from = dedge.first;
    Node to = dedge.second;
    //visit dedges from a leaf as these can be named
    if (G->isLeaf(from)){//from->isLeaf()) {
      if (flat_names.find(dedge) == flat_names.end()) {
        boost::dynamic_bitset<> temp(Gamma_size + 1);
        flat_names[dedge] = temp;
      }
      flat_names.at(dedge)[static_cast<int>( leaf_ids[from->getName()] )] = 1;
      //  flat_names[dedge].insert(leaf_ids[from->getName()]);
      //bl - hack
      long int g_id = set2id(flat_names[dedge]);
      if (G->hasFather(from) and G->getEdgeToFather(from) and G->getEdgeToFather(from)->hasLength()) {//from->hasDistanceToFather())
            Bip_bls[g_id] += G->getEdgeToFather(from)->getLength();
      } else
        Bip_bls[g_id] += 0;
      //mark named
      dedges[dedge] = -1;
      //proceed to dedges in next level -  dedges from cherries can now be named and at least one cherry must exist
      for (auto it_tos = neighbor[to].begin(); it_tos != neighbor[to].end(); it_tos++)
        if ((*it_tos) != from) {
          dedge_type dedge_out;
          dedge_out.first = to;
          dedge_out.second = *it_tos;
          dedges[dedge_out] += 1;
        }
    }
  }

  //dedges identical to usual ale result.

  bool edges_left = false;
  for (auto it = dedges.begin(); it != dedges.end(); it++) {
    dedge_type dedge = (*it).first;
    if (dedges[dedge] != -1)
      edges_left = true;
  }

  if (G->getAllLeaves().size() == 2) {
    Bip_counts[(long int) 1] += weight;

  } else
    while (edges_left) {
      for (auto it = dedges.begin(); it != dedges.end(); it++) {

        dedge_type dedge = (*it).first;
        //Process edges that can be named
        if (dedges[dedge] == 2) {
          Node from = dedge.first;
          Node to = dedge.second;
          vector<dedge_type> dedges_in; //del-loc
          for (auto it_tos = neighbor[from].begin(); it_tos != neighbor[from].end(); it_tos++)
            if (*it_tos != to) {
              dedge_type dedge_in;
              dedge_in.first = *it_tos;
              dedge_in.second = from;
              dedges_in.push_back(dedge_in);
            }
          boost::dynamic_bitset<> leaf_set_in_1 = flat_names[dedges_in[0]];
          boost::dynamic_bitset<> leaf_set_in_2 = flat_names[dedges_in[1]];
          if (flat_names.find(dedge) == flat_names.end()) {
            boost::dynamic_bitset<> temp(Gamma_size + 1);
            flat_names[dedge] = temp;
          }
          flat_names[dedge] = leaf_set_in_1 | leaf_set_in_2;

          long int g_id = set2id(flat_names[dedge]);

          pair<long int, long int> parts;
          long int tmp_id1 = set2id(leaf_set_in_1);
          long int tmp_id2 = set2id(leaf_set_in_2);
          if (tmp_id1 < tmp_id2) {
            parts.first = tmp_id1;
            parts.second = tmp_id2;
          } else {
            parts.first = tmp_id2;
            parts.second = tmp_id1;
          }
          //bl - hack

          if (G->hasFather(from) and G->getFather(from) == to) {
            if (G->getEdgeToFather(from) and G->getEdgeToFather(from)->hasLength())//from->hasDistanceToFather())
              Bip_bls[g_id] += G->getEdgeToFather(from)->getLength();
            else
              Bip_bls[g_id] += 0;
          } else if (G->hasFather(to) and G->getFather(to) == from) {
            if(G->getEdgeToFather(to) and G->getEdgeToFather(to)->hasLength())
                Bip_bls[g_id] += G->getEdgeToFather(to)->getLength();
            else
              Bip_bls[g_id] += 0;
          } else {
            cout << "impossible" << endl;
          }

          //bl - hack
          //INTEGRATED COUNTING
          Dip_counts[g_id][parts] += weight;

          Bip_counts[g_id] += weight;

          //bipartion naming
          if (bip_ids != NULL) bip_ids->insert(g_id);

          //mark named
          dedges[dedge] = -1;
          //proceed to dedges in next level -  new dedges can now be named
          for (auto it_tos = neighbor[to].begin(); it_tos != neighbor[to].end(); it_tos++)
            if ((*it_tos) != from) {
              dedge_type dedge_out;
              dedge_out.first = to;
              dedge_out.second = *it_tos;
              dedges[dedge_out] += 1;
            }
          dedges_in.clear();

        }

      }
      edges_left = false;
      for (auto it = dedges.begin(); it != dedges.end(); it++) {
        dedge_type dedge = (*it).first;
        if (dedges[dedge] != -1) {
          edges_left = true;
          //break;
        }
      }
    }

  //del-locs
  dedges.clear();
  for (auto it = neighbor.begin(); it != neighbor.end(); it++)
    (*it).second.clear();

  neighbor.clear();

  for (auto it = flat_names.begin(); it != flat_names.end(); it++) {
    //(*it).second.clear();
    for (auto i = 0; i < Gamma_size + 1; ++i) { //Resetting all bits
      (*it).second[i] = 0;
      //BipartitionTools::bit0( (*it).second, static_cast<int>(i) );
    }
  }
  flat_names.clear();
  //return return_dips;
}

void approx_posterior::observation(vector<string> trees) {
  scalar_type weight = 1.0;
  for (vector<string>::iterator it = trees.begin(); it != trees.end(); it++) {
    //cout << (*it) << endl;
    decompose(*it, NULL, weight);//del-loc
    observations += weight;
  }

  set_sizes.clear();
  for (map<int, vector<long int> >::iterator it = size_ordered_bips.begin(); it != size_ordered_bips.end(); it++)
    (*it).second.clear();
  size_ordered_bips.clear();

  for (auto it = id_sets.begin(); it != id_sets.end(); it++) {
    size_t size = 0;
    for (auto i = 0; i < Gamma_size + 1; ++i) {
      if ((*it).second[i]) {
        size++;
      }
    }
    set_sizes[(*it).first] = size;
    size_ordered_bips[size].push_back((*it).first);
  }
}



string approx_posterior::mpp_backtrack(long int g_id, map<long int, scalar_type> *qmpp) const {
  //leaf
  if (set_sizes.at(g_id) == 1) {
    stringstream bs;
    bs << Bip_bls.at(g_id) / observations;
    int id = 0;
    for (auto i = 0; i < Gamma_size + 1; ++i) {
      // if ( BipartitionTools::testBit( id_sets.at(g_id), i) ) {
      if (id_sets.at(g_id)[i]) {
        id = i;
        break;
      }
    }

    return id_leaves.at(id) + ":" + bs.str();
  }

  scalar_type max_cp = 0, sum_cp = 0;
  long int max_gp_id = -1;
  long int max_gpp_id = -1;
  for (auto kt = Dip_counts.at(g_id).begin(); kt != Dip_counts.at(g_id).end(); kt++) {
    long int gp_id = (*kt).first.first;
    long int gpp_id = (*kt).first.second;
    scalar_type cp = p_dip(g_id, gp_id, gpp_id) * (*qmpp).at(gp_id) * (*qmpp).at(gpp_id);
    sum_cp += cp;
    if (cp > max_cp) {
      max_cp = cp;
      max_gp_id = gp_id;
      max_gpp_id = gpp_id;
    }
  }
  stringstream bs;
  bs << max_cp / sum_cp << ":" << Bip_bls.at(g_id) / Bip_counts.at(g_id);
  return "(" + mpp_backtrack(max_gp_id, qmpp) + "," + mpp_backtrack(max_gpp_id, qmpp) + ")" + bs.str();
}

vector<string> approx_posterior::all_trees(boost::dynamic_bitset<> gamma) const {
  vector<string> all_trees_g;
  set<int> gamma_s; //not very efficient to use sets again, but this function is rarely used
  for (auto i = 0; i < Gamma_size + 1; ++i) {
    //  if ( BipartitionTools::testBit( gamma, i) ) {
    if (gamma[i]) {
      gamma_s.insert(i);
    }
  }
  if (gamma_s.size() == 1) {
    all_trees_g.push_back(id_leaves.at(*(gamma_s.begin())));
  } else {
    set<set<int> > P_gamma = powerset<set<int> >(gamma_s);//del-loc
    for (set<set<int> >::iterator st = P_gamma.begin(); st != P_gamma.end(); st++) {
      if (gamma_s.size() > (*st).size() and (*st).size() > 0 and (*st).count(*(gamma_s.begin())) == 1) {

        set<int> not_st;//del-loc
        //int* st_bitV = new int[nbint];
        boost::dynamic_bitset<> st_bitV(Gamma_size + 1);
        for (auto it = (*st).begin(); it != (*st).end(); ++it) {
          st_bitV[(*it)] = 1;
        }


        vector<string> all_trees_gp = all_trees(st_bitV);//del-loc

        for (set<int>::iterator nst = gamma_s.begin(); nst != gamma_s.end(); nst++)
          if ((*st).count(*nst) == 0)
            not_st.insert(*nst);

        boost::dynamic_bitset<> not_st_bitV(Gamma_size + 1);
        for (auto it = not_st.begin(); it != not_st.end(); ++it) {
          not_st_bitV[(*it)] = 1;
        }


        vector<string> all_trees_gpp = all_trees(not_st_bitV);//del-loc

        for (vector<string>::iterator lt = all_trees_gp.begin(); lt != all_trees_gp.end(); lt++)
          for (vector<string>::iterator rt = all_trees_gpp.begin(); rt != all_trees_gpp.end(); rt++) {
            all_trees_g.push_back("(" + (*lt) + "," + (*rt) + ")");
          }
        not_st.clear();
        all_trees_gp.clear();
        all_trees_gpp.clear();
      }
    }

    P_gamma.clear();
  }
//  if (Gamma==gamma)
  if (Gamma_size == (int) gamma_s.size())
    for (vector<string>::iterator it = all_trees_g.begin(); it != all_trees_g.end(); it++) {
      (*it) += ";\n";
    }
  return all_trees_g;
}


scalar_type approx_posterior::count_all_trees(boost::dynamic_bitset<> gamma) const {
  scalar_type count_trees_g = 0;
  set<int> gamma_s; //not very efficient to use sets again, but this function is rarely used
  for (auto i = 0; i < Gamma_size + 1; ++i) {
//        if ( BipartitionTools::testBit( gamma, i) ) {
    if (gamma[i]) {
      gamma_s.insert(i);
    }
  }
  if (gamma_s.size() == 1) {
    count_trees_g = 1;
  } else {
    set<set<int> > P_gamma = powerset<set<int> >(gamma_s);//del-loc
    for (set<set<int> >::iterator st = P_gamma.begin(); st != P_gamma.end(); st++) {
      if (gamma_s.size() > (*st).size() and (*st).size() > 0 and (*st).count(*(gamma_s.begin())) == 1) {

        set<int> not_st;//del-loc
        boost::dynamic_bitset<> st_bitV(Gamma_size + 1);
        for (auto it = (*st).begin(); it != (*st).end(); ++it) {
          st_bitV[(*it)] = 1;
        }


        scalar_type count_trees_gp = count_all_trees(st_bitV);//del-loc

        for (set<int>::iterator nst = gamma_s.begin(); nst != gamma_s.end(); nst++)
          if ((*st).count(*nst) == 0)
            not_st.insert(*nst);
        /*
            int* not_st_bitV = new int[nbint];
            for (auto i =0 ; i < nbint  ; ++i) {
                not_st_bitV[i] = 0;
                //BipartitionTools::bit0( not_st_bitV, i) ;
            }
            for (auto it = not_st.begin() ; it != not_st.end()  ; ++it) {
                BipartitionTools::bit1( not_st_bitV, (*it) ) ;
            }
            */
        boost::dynamic_bitset<> not_st_bitV(Gamma_size + 1);
        for (auto it = not_st.begin(); it != not_st.end(); ++it) {
          not_st_bitV[(*it)] = 1;
        }

        scalar_type count_trees_gpp = count_all_trees(not_st_bitV);//del-loc

        count_trees_g += count_trees_gp * count_trees_gpp;
        not_st.clear();
      }
    }

    P_gamma.clear();
  }
  return count_trees_g;
}

scalar_type approx_posterior::count_trees(long int g_id) const {
  scalar_type count_trees_g = 0;
  //std::map <std::set <int>,long int>  set_ids;//del-loc
  //std::map< long int, std::set <int> > id_sets;//del-loc
  //long int g_id=set_ids[gamma];
  boost::dynamic_bitset<> gamma = id_sets.at(g_id);
  int gamma_size = 0;
  for (auto i = 0; i < Gamma_size + 1; ++i) {
    //   if ( BipartitionTools::testBit(gamma, i) ) {
    if (gamma[i]) {
      gamma_size++;
    }
  }
  //gamma.size();
  if (gamma_size == 1) {
    count_trees_g = 1;
  } else {
    set<long int> P_gamma;//=powerset< set<int> >(gamma);//del-loc
    for (auto kt = Dip_counts[g_id].begin(); kt != Dip_counts[g_id].end(); kt++) {
      pair<long int, long int> parts = (*kt).first;
      long int gp_id = parts.first;
      long int gpp_id = parts.second;
      P_gamma.insert(gp_id);
      P_gamma.insert(gpp_id);
    }
    for (auto st = P_gamma.begin(); st != P_gamma.end(); st++) {
      boost::dynamic_bitset<> gammap = id_sets.at((*st));
      int gammap_size = 0;
      bool sameFirstElement = false;
      for (auto i = 0; i < Gamma_size + 1; ++i) {
        //      if ( BipartitionTools::testBit(gammap, i) ) {
        if (gammap[i]) {
          gammap_size++;
        }
      }

      for (auto i = 0; i < Gamma_size + 1; ++i) {
//            if (BipartitionTools::testBit(gamma, i) ) {
        if (gamma[i]) {

//                if ( BipartitionTools::testBit(gammap, i)    ) {
          if (gammap[i]) {

            sameFirstElement = true;
          }
          break;
        }
      }
      if (gamma_size > gammap_size and gammap_size > 0 and sameFirstElement) {
        boost::dynamic_bitset<> not_gammap = ~gammap;
        not_gammap[0] = 0;
        scalar_type count_trees_gp = count_trees((*st));//del-loc
        /*  for (set<int>::iterator nst=gamma.begin();nst!=gamma.end();nst++)
      if (gammap.count(*nst)==0)
        not_gammap.insert(*nst);	*/
        scalar_type count_trees_gpp = count_trees(set_ids.at(not_gammap));//del-loc
        // not_gammap.clear();

        count_trees_g += count_trees_gp * count_trees_gpp;
      }
      // gammap.clear();
    }
    //gamma.clear();
    P_gamma.clear();
  }
  return count_trees_g;
}

