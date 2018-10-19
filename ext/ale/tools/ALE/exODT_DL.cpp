// Modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "exODT_DL.h"

// Bpp includes
#include <Bpp/BppString.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// Treerecs includes
#include <ale/tools/IO/IO.h>
#include <ale/tools/PhyloTreeToolBox.h>

using namespace bpp;
using namespace std;

exODT_DL_model::exODT_DL_model() :
  O_R(1)
{
}




void exODT_DL_model::construct_undated(const string &Sstring) {
  daughter.clear();
  son.clear();
  name_node.clear();
  node_name.clear();
  node_ids.clear();
  S = shared_ptr<tree_type>(IO::newickToPhyloTree(Sstring, true));
  // For each node from the speciestree
  //vector<shared_ptr<bpp::PhyloNode>> nodes = S->getAllNodes();
  auto nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalRecursive_list(*S);

  // Set each branch length to 1.
  for (auto it = nodes.begin(); it != nodes.end(); it++) {
    if (S->hasFather(*it)) {
      auto edge_to_father = S->getEdgeToFather(*it);
      if (edge_to_father) S->getEdgeToFather(*it)->setLength(1);
    }
  }

  // For each node from the speciestree, record node names.
  // name_node gives a node according to the name
  // node_name gives a name accoring to the node.
  // * For leaves, get name.
  // * For internal nodes, gives a node name according to leaves under.
  // \todo node_name is useless until PhyloNode contains a name.
  for (auto it = nodes.begin(); it != nodes.end(); it++) {
    if (S->isLeaf(*it)) {
      name_node[(*it)->getName()] = (*it);
      node_name[(*it)] = (*it)->getName();
    } else {
      vector <string> leafnames = PhyloTreeToolBox::getLeavesNames(*S, (*it));
      sort(leafnames.begin(), leafnames.end());
      stringstream name;
      for (auto leafname = leafnames.begin(); leafname != leafnames.end(); leafname++)
        name << (*leafname) << ".";

      name_node[name.str()] = (*it);
      node_name[(*it)] = name.str();
    }
  }
  daughter.resize(nodes.size());
  son.resize(nodes.size());

  // register species
  last_branch = 0;
  speciesLastLeaf = 0;

  // associate each node to an id (node_ids and i_nodes).
  // this will determines the number of leaves (speciesLastLeaf)
  // and starting to determine the number of branches.
  // \todo Use PhyloTreeToolBox to get direct post-order.
  set<shared_ptr<bpp::PhyloNode>> saw;
  
  for (auto it = name_node.begin(); it != name_node.end(); it++)
    if (S->isLeaf((*it).second)) {
      auto node = (*it).second;
      extant_species[last_branch] = node->getName();
      node_ids[node] = last_branch;
      last_branch++;
      speciesLastLeaf++;
      saw.insert(node);
      // a leaf
      daughter[last_branch] = -1;
      // a leaf
      son[last_branch] = -1;
    }

  //ad-hoc postorder, each internal node is associated to an id and an id to an internal node
  //During the node exploration, set "ID" property to each node with value as the name.
  //
  vector<shared_ptr<bpp::PhyloNode>> next_generation;
  for (auto it = name_node.begin(); it != name_node.end(); it++)
    if (S->isLeaf((*it).second)) {
      auto node = (*it).second;
      next_generation.push_back(node);
    }

  while (next_generation.size()) {
    vector<shared_ptr<bpp::PhyloNode>> new_generation;
    for (auto it = next_generation.begin(); it != next_generation.end(); it++) {
      auto node = (*it);
      if (S->hasFather(node)) {
        auto father = S->getFather(node);
        auto sons = S->getSons(father);
        decltype(node) sister;

        if (sons[0] == node) 
          sister = sons[1];
        else 
          sister = sons[0];

        if (not node_ids.count(father) and saw.count(sister)) {
          node_ids[father] = last_branch;
          last_branch++;
          saw.insert(father);
          new_generation.push_back(father);
        }
      }
    }
    next_generation.clear();
    for (auto it = new_generation.begin();
         it != new_generation.end(); it++)
      next_generation.push_back((*it));
  }

  // Fill daughter and son maps. Daughter contains son left and son son right of a node.
  for (auto it = name_node.begin(); it != name_node.end(); it++) {
    if (not S->isLeaf(it->second))//not (*it).second->isLeaf())
    {
      auto node = (*it).second;
      auto sons = S->getSons(node);//node->getSons();
      daughter[node_ids[node]] = node_ids[sons[0]];
      son[node_ids[node]] = node_ids[sons[1]];
    }
  }
  last_rank = last_branch;
}

void exODT_DL_model::calculate_undatedEs() {
  uE = vector<scalar_type>(last_branch, 0.0);
  for (int e = 0; e < speciesLastLeaf; ++e) {
    scalar_type a = PD;
    scalar_type b = -1.0;
    scalar_type c = PL;
    uE[e] = (-b - sqrt(b * b - 4 * a * c)) / (2.0 * a);
  }
  for (int e = speciesLastLeaf; e < last_branch; ++e) {
    int f = daughter[e];
    int g = son[e];
    scalar_type a = PD;
    scalar_type b = -1.0;
    scalar_type c = PL + PS * uE[f]  * uE[g];
    uE[e] = (-b - sqrt(b * b - 4 * a * c)) / (2.0 * a);
  }
}

void  exODT_DL_model::step_one(shared_ptr<approx_posterior> ale) {
  for (auto it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (auto jt = (*it).second.begin(); jt != (*it).second.end(); jt++) {
      g_ids.push_back((*jt));
      g_id_sizes.push_back((*it).first);
    }
  //root bipartition needs to be handled separately
  g_ids.push_back(-1);
  g_id_sizes.push_back(ale->Gamma_size);
}

void exODT_DL_model::gene_species_mapping(shared_ptr<approx_posterior> ale)
{
  // gene<->species mapping
  if (gid_sps.size() == 0) // If the mapping has not been done yet
  {
    //Test that the species associated to genes are really in the species tree
    set<string> species_set;
    for (auto iter = extant_species.begin(); iter != extant_species.end(); ++iter) {
      species_set.insert(iter->second);
    }

    for (int i = 0; i < (int) g_ids.size(); i++) {
      long int g_id = g_ids[i];

      if (g_id_sizes[i] == 1) {
        int id = 0;
        for (auto i = 0; i < ale->Gamma_size + 1; ++i) {
          if (ale->id_sets[g_id][i]) {
            id = i;
            break;
          }
        }

        string gene_name = ale->id_leaves[id];
        // Set associated species in gid_sps[g_id]
        string species_name = speciesGeneMap.getAssociatedSpecies(gene_name);
        gid_sps[g_id] = species_name;
        if (species_set.find(species_name) == species_set.end()) {
          cout << "Error: gene name " << gene_name << " is associated to species name " << species_name
                    << " that cannot be found in the species tree." << endl;
          exit(-1);
        }
      }
    }
  }
}

void exODT_DL_model::inner_loop(bool g_is_a_leaf,
                             long int &g_id,
                             vector<int> &gp_is,
                             vector<long int> &gpp_is,
                             int i) {
  for (int e = 0; e < last_branch; e++) {
    bool s_is_leaf = (e < speciesLastLeaf);
    int f = 0;
    int g = 0;
    if (not s_is_leaf) {
      f = daughter[e];
      g = son[e];
    }
    scalar_type uq_sum = 0;
    // S leaf and G leaf
    if (e < speciesLastLeaf and g_is_a_leaf and extant_species[e] == gid_sps[g_id]) {
      // present
      uq_sum += PS;
    }
    // G internal
    if (not g_is_a_leaf) {
      int N_parts = gp_is.size();
      for (int i = 0; i < N_parts; i++) {
        int gp_i = gp_is[i];
        int gpp_i = gpp_is[i];
        if (not s_is_leaf) {
          uq_sum += PS * (uq[gp_i][f] * uq[gpp_i][g] + uq[gp_i][g] * uq[gpp_i][f]);
        }
        // D event
        uq_sum += PD * (uq[gp_i][e] * uq[gpp_i][e] * 2);
      }
    }
    if (not s_is_leaf) {
      // SL event
      uq_sum += PS * (uq[i][f] * uE[g] + uq[i][g] * uE[f]);
    }
    uq[i][e] = uq_sum / (1.0 - 2.0 * PD * uE[e]);
  }
}

scalar_type exODT_DL_model::pun(shared_ptr<approx_posterior> ale, bool verbose) {
  scalar_type survive = 0;
  scalar_type root_sum = 0;
  scalar_type O_norm = 0;
  
  g_ids.clear();
  g_id_sizes.clear();
  step_one(ale);
  gene_species_mapping(ale);

  for (int i = 0; i < (int) g_ids.size(); i++) {
    long int g_id = g_ids[i];
    g_id2i[g_id] = i;
  }
  vector<scalar_type> zeros(last_branch, 0.0);
  uq = vector<vector<scalar_type>>(g_ids.size(),zeros);

  for (int i = 0; i < (int) g_ids.size(); i++) {
    // directed partition (dip) gamma's id
    long int g_id = g_ids[i];
    bool is_a_leaf = (g_id_sizes[i] == 1);
    vector<int> gp_is;
    vector<long int> gpp_is;
    if (g_id != -1) {
      for (unordered_map<pair<long int, long int>, scalar_type>::iterator kt = ale->Dip_counts[g_id].begin();
          kt != ale->Dip_counts[g_id].end(); kt++) {
        pair<long int, long int> parts = (*kt).first;
        long int gp_id = parts.first;
        long int gpp_id = parts.second;
        int gp_i = g_id2i[parts.first];
        int gpp_i = g_id2i[parts.second];
        gp_is.push_back(gp_i);
        gpp_is.push_back(gpp_i);
      }
    } else {
      //XX
      //root bipartition needs to be handled separately
      map<set<long int>, int> bip_parts;
      for (map<long int, scalar_type>::iterator it = ale->Bip_counts.begin();
          it != ale->Bip_counts.end(); it++) {
        long int gp_id = (*it).first;
        boost::dynamic_bitset<> gamma = ale->id_sets.at(gp_id);
        boost::dynamic_bitset<> not_gamma = ~gamma;
        not_gamma[0] = 0;
        long int gpp_id = ale->set_ids.at(not_gamma);

        set<long int> parts;
        parts.insert(gp_id);
        parts.insert(gpp_id);
        bip_parts[parts] = 1;
      }
      for (map<set<long int>, int>::iterator kt = bip_parts.begin(); kt != bip_parts.end(); kt++) {
        vector<long int> parts;
        for (set<long int>::iterator sit = (*kt).first.begin(); sit != (*kt).first.end(); sit++) {
          parts.push_back((*sit));
        }
        long int gp_id = parts[0];
        int gp_i = g_id2i[parts[0]];
        int gpp_i = g_id2i[parts[1]];
        gp_is.push_back(gp_i);
        gpp_is.push_back(gpp_i);
      }
      bip_parts.clear();
    }
    inner_loop(is_a_leaf, g_id, gp_is, gpp_is, i);
  }
  survive = 0;
  root_sum = 0;
  O_norm = 0;
  for (int e = 0; e < last_branch; e++) {
    scalar_type O_p = 1;
    if (e == (last_branch - 1)) 
      O_p = O_R;
    O_norm += O_p;
    root_sum += uq[g_ids.size() - 1][e] * O_p;
    survive += (1 - uE[e]);
  }
  return root_sum / survive / O_norm * (last_branch);
}
  
void exODT_DL_model::setRates(scalar_type dupRate, scalar_type lossRate) {
    PD = dupRate;
    PL = lossRate;
    PS = 1;
    scalar_type sum = PD + PL + PS;
    PD /= sum;
    PL /= sum;
    PS /= sum;
  }

exODT_DL_model::~exODT_DL_model() { }

