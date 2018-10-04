// Modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "exODT.h"
#include "fractionMissing.h"

// Bpp includes
#include <Bpp/BppString.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// Treerecs includes
#include <ale/tools/IO/IO.h>
#include <ale/tools/PhyloTreeToolBox.h>

using namespace bpp;
using namespace std;

exODT_model::exODT_model() {
  //some default parameters
  string_parameter["BOOTSTRAP_LABELS"] = "no";
  string_parameter["gene_name_separators"] = "_@";
  scalar_parameter["species_field"] = 0;
  scalar_parameter["event_node"] = 0;
  scalar_parameter["min_bip_count"] = -1;
  scalar_parameter["min_branch_lenghts"] = 0;
  // length of "stem" branch above root
  scalar_parameter["stem_length"] = 1;
  //number of discretization slices (subslices) per time slice
  scalar_parameter["D"] = 3;
  scalar_parameter["grid_delta_t"] = 0.005;
  scalar_parameter["min_D"] = 3;
  //number of subdiscretizations for ODE calculations
  scalar_parameter["DD"] = 10;

}


void exODT_model::set_model_parameter(std::string name, std::string value) {
  string_parameter[name] = value;
}


void exODT_model::set_model_parameter(std::string name, scalar_type value) {

  if (name == "delta" or name == "tau" or name == "lambda") {
    if (name == "tau") {
      transfers = (value != 0.0);
    }
    scalar_type N = vector_parameter["N"][0];
    vector_parameter[name].clear();
    for (int branch = 0; branch < last_branch; branch++)
      if (name == "tau")
        vector_parameter[name].push_back(value / N);
      else
        vector_parameter[name].push_back(value);
    if (name == "tau")
      scalar_parameter[name + "_avg"] = value / N;
    else
      scalar_parameter[name + "_avg"] = value;

  } else if (name == "N" or name == "Delta_bar" or name == "Lambda_bar") {
    vector_parameter[name].clear();
    for (int rank = 0; rank < last_rank; rank++)
      vector_parameter[name].push_back(value);
  } else
    scalar_parameter[name] = value;
}


static scalar_type EPSILON = std::numeric_limits<scalar_type>::min();

void exODT_model::construct_undated(const std::string &Sstring, const std::string &fractionMissingFile) {
  daughter.clear();
  son.clear();
  name_node.clear();
  node_name.clear();
  node_ids.clear();
  id_nodes.clear();

  string_parameter["S_un"] = Sstring;
  S = std::shared_ptr<tree_type>(IO::newickToPhyloTree(string_parameter["S_un"], true));
  // For each node from the speciestree
  //std::vector<std::shared_ptr<bpp::PhyloNode>> nodes = S->getAllNodes();
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
      std::vector <std::string> leafnames = PhyloTreeToolBox::getLeavesNames(*S, (*it));
      std::sort(leafnames.begin(), leafnames.end());
      std::stringstream name;
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
  last_leaf = 0;

  // associate each node to an id (node_ids and i_nodes).
  // this will determines the number of leaves (last_leaf)
  // and starting to determine the number of branches.
  // \todo Use PhyloTreeToolBox to get direct post-order.
  std::set<std::shared_ptr<bpp::PhyloNode>> saw;
  
  for (auto it = name_node.begin(); it != name_node.end(); it++)
    if (S->isLeaf((*it).second)) {
      auto node = (*it).second;
      extant_species[last_branch] = node->getName();
      node_ids[node] = last_branch;
      id_nodes[last_branch] = node;
      last_branch++;
      last_leaf++;
      saw.insert(node);
      // a leaf
      daughter[last_branch] = -1;
      // a leaf
      son[last_branch] = -1;
    }

  //ad-hoc postorder, each internal node is associated to an id and an id to an internal node
  //During the node exploration, set "ID" property to each node with value as the name.
  //
  std::vector<std::shared_ptr<bpp::PhyloNode>> next_generation;
  for (auto it = name_node.begin(); it != name_node.end(); it++)
    if (S->isLeaf((*it).second)) {
      auto node = (*it).second;
      next_generation.push_back(node);
    }

  while (next_generation.size()) {
    std::vector<std::shared_ptr<bpp::PhyloNode>> new_generation;
    for (auto it = next_generation.begin(); it != next_generation.end(); it++) {
      auto node = (*it);
      if (S->hasFather(node)) {
        auto father = S->getFather(node);
        auto sons = S->getSons(father);
        decltype(node) sister;

        if (sons[0] == node) sister = sons[1];
        else sister = sons[0];

        if (not node_ids.count(father) and saw.count(sister)) {
          node_ids[father] = last_branch;
          id_nodes[last_branch] = father;
          std::stringstream name;
          name << last_branch;
          if (S->hasFather(father) and S->getEdgeToFather(father))
              S->getEdgeToFather(father)->setProperty("ID", BppString(name.str()));

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
  // End of the ad-hoc post-order constitution

  // Collect each bootstrap value (from upper branch of a given node).
  // If the edge has no bootstrap, set default value as -1.
  // This, fills rank2label where rank is the branch and label the bootstrap.
  // Each branch is edited with a property "ID".
  // \todo Each label of rank2label, a bootstrap is an integer. Sometimes we can find bootstrap as floating values.
  for (auto it = node_ids.begin(); it != node_ids.end(); it++) {
    auto node = (*it).first;
    int branch = (*it).second;
    std::stringstream out;
    std::stringstream out1;
    std::stringstream out2;
    out1 << t_begin[branch];
    out2 << t_end[branch];
    int rank = branch;
    out << rank;
    if (S->hasFather(node) and S->getEdgeToFather(node) and (S->getEdgeToFather(node)->hasBootstrapValue())) {
      rank2label[rank] = S->getEdgeToFather(node)->getBootstrapValue();
    } else {
      rank2label[rank] = -1;
    }
    if (S->hasFather(node) and S->getEdgeToFather(node))
      S->getEdgeToFather(node)->setProperty("ID", BppString(out.str()));
  }

  // Save the resulting speciestree with the "ID" label at each branch.
  string_parameter["S_with_ranks"] = IO::PhyloTreeToNewick(*S, false, "ID");

  //map <string,map<string,int> > ancestral_names;
  //map <int,map<int,int> > ancestral;

  // Initialize the ancestor matrix (nb_branches x bn_branches) with zeros.
  // Value take one if node has parent's relationship with an other.
  ancestors.clear();
  for (int e = 0; e < last_branch; e++) {
    std::vector<int> tmp;
    ancestors.push_back(tmp);
    for (int f = 0; f < last_branch; f++)
      ancestral[e][f] = 0;
  }

  // For each node in the species tree, fill ancestral_names and ancestral which are both binary matrices.
  for (auto it = nodes.begin(); it != nodes.end(); it++) {
    auto node = (*it);
    int edge = node_ids[node];
    std::stringstream name_from; //contains node name or edge if it's root.
    if (edge < last_leaf)
      name_from << node_name[node];
    else
      name_from << edge;

    std::map<std::string, int> tmp;
    ancestral_names[name_from.str()] = tmp;
    while (S->hasFather(node)) { //While we are not at the root.
      std::stringstream name_to;
      int f = node_ids[node];
      if (f < last_leaf)
        name_to << node_name[node];
      else
        name_to << f;

      node = S->getFather(node);
      ancestral_names[name_from.str()][name_to.str()] = 1;
      if (not ancestral[edge][f])
        ancestors[edge].push_back(f);

      ancestral[edge][f] = 1;
    }
    std::stringstream name_to;
    int f = node_ids[node];
    name_to << f;
    ancestral_names[name_from.str()][name_to.str()] = 1;
    if (not ancestral[edge][f])
      ancestors[edge].push_back(f);
    ancestral[edge][f] = 1;
  }

  // Fill daughter and son maps. Daughter contains son left and son son right of a node.
  for (auto it = name_node.begin(); it != name_node.end(); it++) {
    if (not S->isLeaf(it->second))//not (*it).second->isLeaf())
    {
      auto node = (*it).second;
      auto sons = S->getSons(node);//node->getSons();
      daughter[node_ids[node]] = node_ids[sons[0]];
      son[node_ids[node]] = node_ids[sons[1]];
      //cout << node_ids[node] << " => " << node_ids[sons[0]] << " & " << node_ids[sons[1]] << endl;
      //cout << node_name[node] << " => " << node_name[sons[0]] << " & " << node_name[sons[1]] << endl;

    }
  }

  //Init branch_counts
  branch_counts["Os"].clear();
  branch_counts["Ds"].clear();
  branch_counts["Ts"].clear();
  branch_counts["Tfroms"].clear();
  branch_counts["Ls"].clear();
  branch_counts["count"].clear();
  branch_counts["copies"].clear();
  branch_counts["singleton"].clear();

  for (int e = 0; e < last_branch; e++) {
    branch_counts["Os"].push_back(0);
    branch_counts["Ds"].push_back(0);
    branch_counts["Ts"].push_back(0);
    branch_counts["Tfroms"].push_back(0);
    branch_counts["Ls"].push_back(0);
    branch_counts["count"].push_back(0);
    branch_counts["copies"].push_back(0);
    branch_counts["singleton"].push_back(0);

  }
  T_to_from.clear();
  // For each edge, init T_to_from
  for (int e = 0; e < last_branch; e++) {
    std::vector<scalar_type> tmp;
    T_to_from.push_back(tmp);
    for (int f = 0; f < last_branch; f++)
      T_to_from[e].push_back(0);
  }


  last_rank = last_branch;
  set_model_parameter("N", 1);


  //Put default values for the fraction of missing genes at the leaves.
  vector_parameter["fraction_missing"] = std::vector<scalar_type>(last_leaf, 0.0);
  //Put user-defined values, if available
  if (fractionMissingFile == "") {

  } else {
    fraction_missing = readFractionMissingFile(fractionMissingFile);
    // Now we need to fill up the vector_parameter, and we have to be careful about the order.
    std::size_t index = 0;
    for (auto it = name_node.begin(); it != name_node.end(); it++) {
      if (S->isLeaf(it->second))//(*it).second->isLeaf())
      {
        auto node = (*it).second;
        std::string currentSpecies = node->getName();
        vector_parameter["fraction_missing"][index] = fraction_missing[currentSpecies];
        index++;
      }
    }
    VectorTools::print(vector_parameter["fraction_missing"]);
  }
}

void exODT_model::calculate_undatedEs() {
  uE.clear();
  fm.clear();
  mPTE_ancestral_correction.clear();
  PD.clear();
  PT.clear();
  PL.clear();
  PS.clear();
  scalar_type P_T = 0;
  for (int f = 0; f < last_branch; f++)
    P_T += vector_parameter["tau"][f] / (scalar_type) last_branch;

  for (int e = 0; e < last_branch; e++) {
    //cout << e << " is " << node_name[id_nodes[e]] << endl;
    scalar_type P_D = vector_parameter["delta"][e];
    scalar_type P_L = vector_parameter["lambda"][e];
    scalar_type P_S = 1;
    scalar_type tmp = P_D + P_T + P_L + P_S;
    P_D /= tmp;
    P_L /= tmp;
    P_S /= tmp;
    PD.push_back(P_D);
    PT.push_back(P_T / tmp);
    PL.push_back(P_L);
    PS.push_back(P_S);
    uE.push_back(0);
    if (e < last_leaf) { // we are at a leaf
      fm.push_back(vector_parameter["fraction_missing"][e]);
    } else {
      fm.push_back(0);
    }
    mPTE_ancestral_correction.push_back(0);
  }


// In the loop below with 4 iterations, we calculate the mean probability mPTE for a gene to become extinct across all branches.
  mPTE = 0;
  for (int i = 0; i < 4; i++) {
    scalar_type newmPTE = 0;
    //vector<scalar_type> ancestral_correction;
    if (i > 0) // There should be no need for this loop at the first iteration, because then it leaves mPTE_ancestral_correction at 0.
    {
      for (int edge = 0; edge < last_branch; edge++) {
        mPTE_ancestral_correction[edge] = 0;
        //for (map<int,int>::iterator it=ancestral[edge].begin();( it!=ancestral[edge].end() and i>0);it++)
        for (auto it = ancestors[edge].begin(); it != ancestors[edge].end(); it++) {
          //int f=(*it).first;
          int f = (*it);
          //if (ancestral[edge][f]==1)
          mPTE_ancestral_correction[edge] +=
              (PT[f] / (scalar_type) last_branch) * uE[f]; //That's how we forbid transfers to ancestors of a branch
        }
      }
    }

    for (int e = 0; e < last_branch; e++) {
      if (e < last_leaf) // we are at a leaf, there cannot be a speciation event, but the gene has been lost
      {
        uE[e] = PL[e] + PD[e] * uE[e] * uE[e] + uE[e] * (mPTE - mPTE_ancestral_correction[e]);
      } else // Not at a leaf: the gene was lost once on branch e, or on all descendants, including after speciation
      {
        int f = daughter[e];
        int g = son[e];
        uE[e] = PL[e] + PS[e] * uE[f] * uE[g] + PD[e] * uE[e] * uE[e] + uE[e] * (mPTE - mPTE_ancestral_correction[e]);
      }
      newmPTE += (PT[e] / (scalar_type) last_branch) * uE[e];
    }
    mPTE = newmPTE;
  } // End of the loop to compute mPTE

  // Now we add one more update of mPTE to take into account the fraction of missing genes, which had been ignored so far.
  scalar_type newmPTE = 0;
  for (int e = 0; e < last_branch; e++) {
    if (e < last_leaf) // we are at a leaf: either the gene has been lost, or it is one of those missing genes
    {
      uE[e] = (1 - fm[e]) * uE[e] + fm[e];
    } else // Not at a leaf: the gene was lost once on branch e, or on all descendants
    {
      int f = daughter[e];
      int g = son[e];
      uE[e] = PL[e] + PS[e] * uE[f] * uE[g] + PD[e] * uE[e] * uE[e] + uE[e] * (mPTE - mPTE_ancestral_correction[e]);
    }
    newmPTE += (PT[e] / (scalar_type) last_branch) * uE[e];
  }
  mPTE = newmPTE;

}

void  exODT_model::clear_all()
{
  mPTuq_ancestral_correction.clear();
  uq.clear();
  mPTuq.clear();//XX
  for (auto it = q.begin(); it != q.end(); it++) {
    for (auto jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      (*jt).second.clear();
    (*it).second.clear();
  }
  q.clear();

  //directed partitions and their sizes
  //vector <long int>  g_ids;
  //vector <long int>  g_id_sizes;
  g_ids.clear();
  g_id_sizes.clear();
}

void  exODT_model::step_one(std::shared_ptr<approx_posterior> ale) {
  for (auto it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (auto jt = (*it).second.begin(); jt != (*it).second.end(); jt++) {
      g_ids.push_back((*jt));
      g_id_sizes.push_back((*it).first);
    }
  //root bipartition needs to be handled separately
  g_ids.push_back(-1);
  g_id_sizes.push_back(ale->Gamma_size);

  root_i = g_ids.size() - 1;
}

void exODT_model::gene_secies_mapping(std::shared_ptr<approx_posterior> ale)
{

  // gene<->species mapping
  if (gid_sps.size() == 0) // If the mapping has not been done yet
  {
    //Test that the species associated to genes are really in the species tree
    std::set<std::string> species_set;
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

        std::string gene_name = ale->id_leaves[id];

//        std::vector<std::string> tokens = Utils::splitString(gene_name, string_parameter["gene_name_separators"].c_str());
//        std::string species_name;
//        if ((int) scalar_parameter["species_field"] == -1) {
//          species_name = tokens[1];
//          for (int fi = 2; fi < tokens.size(); fi++)
//            species_name += "_" + tokens[fi];
//        } else
//          species_name = tokens[(int) scalar_parameter["species_field"]];


        // Set associated species in gid_sps[g_id]
        std::string species_name = speciesGeneMap.getAssociatedSpecies(gene_name);
        gid_sps[g_id] = species_name;
        if (species_set.find(species_name) == species_set.end()) {
          std::cout << "Error: gene name " << gene_name << " is associated to species name " << species_name
                    << " that cannot be found in the species tree." << std::endl;
          exit(-1);
        }
      }
    }
  }
}

void exODT_model::inner_loop(std::shared_ptr<approx_posterior> ale, bool is_a_leaf,
                             long int &g_id,
                             std::vector<int> &gp_is,
                             std::vector<long int> &gpp_is,
                             std::vector<scalar_type> &p_part,
                             int i,
                             scalar_type &new_mPTuq) {
  for (int e = 0; e < last_branch; e++) {
    int f = 0;
    int g = 0;
    if (not(e < last_leaf)) {
      f = daughter[e];
      g = son[e];
    }
    scalar_type uq_sum = 0;
    // S leaf and G leaf
    if (e < last_leaf and is_a_leaf and extant_species[e] == gid_sps[g_id]) {
      // present
      uq_sum += PS[e] * 1;
    }
    // G internal
    if (not is_a_leaf) {
      int N_parts = gp_is.size();
      for (int i = 0; i < N_parts; i++) {
        int gp_i = gp_is[i];
        int gpp_i = gpp_is[i];
        scalar_type pp = p_part[i];
        if (not(e < last_leaf)) {
          uq_sum += PS[e] * (uq[gp_i][f] * uq[gpp_i][g] + uq[gp_i][g] * uq[gpp_i][f]) * pp;
        }
        // D event
        uq_sum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2) * pp;
        // T event
        if (transfers) {
          uq_sum += (uq[gp_i][e] * (mPTuq[gpp_i] - mPTuq_ancestral_correction[gpp_i][e]) +
                    uq[gpp_i][e] * (mPTuq[gp_i] - mPTuq_ancestral_correction[gp_i][e])) * pp;
        }
      }
    }
    if (not(e < last_leaf)) {
      // SL event
      uq_sum += PS[e] * (uq[i][f] * uE[g] + uq[i][g] * uE[f]);
    }
    // DL event
    uq_sum += PD[e] * (uq[i][e] * uE[e] * 2);
    // TL event
    if (transfers) {
      uq_sum += ((mPTuq[i] - mPTuq_ancestral_correction[i][e]) * uE[e] +
               uq[i][e] * (mPTE - mPTE_ancestral_correction[e]));
    }
    if (uq_sum < EPSILON) uq_sum = EPSILON;
    uq[i][e] = uq_sum;
    if (transfers) {
      new_mPTuq += (PT[e] / (scalar_type) last_branch) * uq_sum;
      mPTuq_ancestral_correction[i][e] = 0;
      for (std::vector<int>::iterator it = ancestors[e].begin(); it != ancestors[e].end(); it++) {
        mPTuq_ancestral_correction[i][e] += (PT[*it] / (scalar_type) last_branch) * uq_sum;
      }
    }
  }
  mPTuq[i] = new_mPTuq;

}

scalar_type exODT_model::pun(std::shared_ptr<approx_posterior> ale, bool verbose) {
  scalar_type survive = 0;
  scalar_type root_sum = 0;
  scalar_type O_norm = 0;

  ale_pointer = ale;
  clear_all();
  step_one(ale);
  gene_secies_mapping(ale);


  //map <long int, int> g_id2i;
  //XX ancestral_correction ..
  for (int i = 0; i < (int) g_ids.size(); i++) {
    long int g_id = g_ids[i];
    g_id2i[g_id] = i;

    if (not(i < (int) uq.size())) {
      std::vector<scalar_type> tmp;
      uq.push_back(tmp);
      std::vector<scalar_type> tmp2;
      mPTuq_ancestral_correction.push_back(tmp2);

      mPTuq.push_back(0);
    } else
      mPTuq[i] = 0;

    for (int e = 0; e < last_branch; e++)
      if (not(e < (int) uq[i].size())) {
        uq[i].push_back(0);
        mPTuq_ancestral_correction[i].push_back(0);
      } else {
        uq[i][e] = 0;
        mPTuq_ancestral_correction[i][e] = 0;
      }
  }

  for (int iter = 0; iter < 4; iter++) {
    for (int i = 0; i < (int) g_ids.size(); i++) {

      scalar_type new_mPTuq = 0;

      // directed partition (dip) gamma's id
      bool is_a_leaf = false;
      long int g_id = g_ids[i];
      if (g_id_sizes[i] == 1)
        is_a_leaf = true;
      std::vector<int> gp_is;
      std::vector<long int> gpp_is;
      std::vector<scalar_type> p_part;
      if (g_id != -1)
        for (std::unordered_map<std::pair<long int, long int>, scalar_type>::iterator kt = ale->Dip_counts[g_id].begin();
             kt != ale->Dip_counts[g_id].end(); kt++) {
          std::pair<long int, long int> parts = (*kt).first;
          long int gp_id = parts.first;
          long int gpp_id = parts.second;
          int gp_i = g_id2i[parts.first];
          int gpp_i = g_id2i[parts.second];
          gp_is.push_back(gp_i);
          gpp_is.push_back(gpp_i);
          if (ale->Bip_counts[g_id] <= scalar_parameter["min_bip_count"])
            p_part.push_back(0);
          else
            p_part.push_back(
                pow((scalar_type) ale->p_dip(g_id, gp_id, gpp_id), (scalar_type) scalar_parameter["seq_beta"]));//set pp
        }
      else {
        //XX
        //root bipartition needs to be handled separately
        std::map<std::set<long int>, int> bip_parts;
        for (std::map<long int, scalar_type>::iterator it = ale->Bip_counts.begin();
             it != ale->Bip_counts.end(); it++) {
          long int gp_id = (*it).first;
          boost::dynamic_bitset<> gamma = ale->id_sets.at(gp_id);
          boost::dynamic_bitset<> not_gamma = ~gamma;
          not_gamma[0] = 0;
          long int gpp_id = ale->set_ids.at(not_gamma);

          std::set<long int> parts;
          parts.insert(gp_id);
          parts.insert(gpp_id);
          bip_parts[parts] = 1;
          // gamma.clear();
          // not_gamma.clear();
        }
        for (std::map<std::set<long int>, int>::iterator kt = bip_parts.begin(); kt != bip_parts.end(); kt++) {
          std::vector<long int> parts;
          for (std::set<long int>::iterator sit = (*kt).first.begin(); sit != (*kt).first.end(); sit++) {
            parts.push_back((*sit));
          }
          long int gp_id = parts[0];
          //long int gpp_id=parts[1];

          int gp_i = g_id2i[parts[0]];
          int gpp_i = g_id2i[parts[1]];
          gp_is.push_back(gp_i);
          gpp_is.push_back(gpp_i);

          //Here we can create a new ale->Bip_counts[gp_id], in particular for leaves.
          //We may want to add the leaf entries for Bip_counts when Bip_counts is first created.
          if ((ale->Bip_counts[gp_id] <= scalar_parameter.at("min_bip_count")) and not (ale->Gamma_size < 4))
            p_part.push_back(0);
          else
            p_part.push_back(pow((scalar_type) ale->p_bip(gp_id), (scalar_type) scalar_parameter["seq_beta"]));//set pp
        }
        bip_parts.clear();
      }
      //######################################################################################################################
      //#########################################INNNER LOOP##################################################################
      //######################################################################################################################

      inner_loop(ale, is_a_leaf, g_id, gp_is, gpp_is, p_part, i, new_mPTuq);
      //######################################################################################################################
      //#########################################INNNER LOOP##################################################################
      //######################################################################################################################
    }
    survive = 0;
    root_sum = 0;
    O_norm = 0;
    for (int e = 0; e < last_branch; e++) {
      scalar_type O_p = 1;
      if (e == (last_branch - 1)) O_p = scalar_parameter["O_R"];
      O_norm += O_p;
      root_sum += uq[root_i][e] * O_p;
      survive += (1 - uE[e]);
    }
    //cout << root_sum/survive << endl;
  }
  return root_sum / survive / O_norm * (last_branch);

}

std::string exODT_model::sample_undated() {

  scalar_type r = bpp::RandomTools::giveRandomNumberBetweenZeroAndEntry(1);

  scalar_type root_sum = 0;

  for (int e = 0; e < last_branch; e++) {
    scalar_type O_p = 1;
    if (e == (last_branch - 1)) O_p = scalar_parameter["O_R"];
    root_sum += uq[g_ids.size() - 1][e] * O_p + EPSILON;
  }
  scalar_type root_resum = 0;
  //scalar_type O_norm = 0;

  for (int e = 0; e < last_branch; e++) {
    scalar_type O_p = 1;
    if (e == (last_branch - 1)) O_p = scalar_parameter["O_R"];

    root_resum += uq[root_i][e] * O_p + EPSILON;
    if (r * root_sum < root_resum) {
      register_O(e);
      return sample_undated(e, root_i, "O") + ";";
    }
  }
  return "-!=-";
}

std::string exODT_model::sample_undated(int e, int i, std::string last_event, std::string branch_string) {


  scalar_type r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);

  bool is_a_leaf = false;
  long int g_id = g_ids[i];
  if (g_id_sizes[i] == 1)
    is_a_leaf = true;

  std::stringstream bl;
  if (ale_pointer->Bip_counts.count(g_id) and ale_pointer->Bip_counts[g_id] > 0)
    bl << std::max(ale_pointer->Bip_bls[g_id] / ale_pointer->Bip_counts[g_id],
                   (scalar_type) scalar_parameter["min_branch_lenghts"]);
  else
    bl << std::max(ale_pointer->Bip_bls[g_id] / ale_pointer->observations,
                   (scalar_type) scalar_parameter["min_branch_lenghts"]);
  std::string branch_length = bl.str();

  std::vector<int> gp_is;
  std::vector<long int> gpp_is;
  std::vector<scalar_type> p_part;
  if (g_id != -1)
    for (std::unordered_map<std::pair<long int, long int>, scalar_type>::iterator kt = ale_pointer->Dip_counts[g_id].begin();
         kt != ale_pointer->Dip_counts[g_id].end(); kt++) {
      std::pair<long int, long int> parts = (*kt).first;
      long int gp_id = parts.first;
      long int gpp_id = parts.second;
      int gp_i = g_id2i[parts.first];
      int gpp_i = g_id2i[parts.second];
      gp_is.push_back(gp_i);
      gpp_is.push_back(gpp_i);
      if (ale_pointer->Bip_counts[g_id] <= scalar_parameter["min_bip_count"])
        p_part.push_back(0);
      else
        p_part.push_back(pow((scalar_type) ale_pointer->p_dip(g_id, gp_id, gpp_id),
                             (scalar_type) scalar_parameter["seq_beta"]));//set pp
    }
  else {
    //root bipartition needs to be handled separately
    std::map<std::set<long int>, int> bip_parts;
    for (std::map<long int, scalar_type>::iterator it = ale_pointer->Bip_counts.begin();
         it != ale_pointer->Bip_counts.end(); it++) {
      long int gp_id = (*it).first;
      boost::dynamic_bitset<> gamma = ale_pointer->id_sets.at(gp_id);
      boost::dynamic_bitset<> not_gamma = ~gamma;
      not_gamma[0] = 0;
      long int gpp_id = ale_pointer->set_ids.at(not_gamma);

      std::set<long int> parts;
      parts.insert(gp_id);
      parts.insert(gpp_id);
      bip_parts[parts] = 1;
      // gamma.clear();
      // not_gamma.clear();
    }
    for (std::map<std::set<long int>, int>::iterator kt = bip_parts.begin(); kt != bip_parts.end(); kt++) {
      std::vector<long int> parts;
      for (std::set<long int>::iterator sit = (*kt).first.begin(); sit != (*kt).first.end(); sit++) {
        parts.push_back((*sit));
      }
      long int gp_id = parts[0];
      //long int gpp_id=parts[1];

      int gp_i = g_id2i[parts[0]];
      int gpp_i = g_id2i[parts[1]];
      gp_is.push_back(gp_i);
      gpp_is.push_back(gpp_i);

      //Here we can create a new ale->Bip_counts[gp_id], in particular for leaves.
      //We may want to add the leaf entries for Bip_counts when Bip_counts is first created.
      if ((ale_pointer->Bip_counts[gp_id] <= scalar_parameter.at("min_bip_count")) and not (ale_pointer->Gamma_size < 4))
        p_part.push_back(0);
      else
        p_part.push_back(
            pow((scalar_type) ale_pointer->p_bip(gp_id), (scalar_type) scalar_parameter["seq_beta"]));//set pp
    }
    bip_parts.clear();
  }


  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################
  scalar_type uq_sum = 0;
  // S leaf and G leaf
  if (e < last_leaf and is_a_leaf and extant_species[e] == gid_sps[g_id]) {
    // present
    uq_sum += PS[e] * 1 + EPSILON;
  }

  // G internal
  if (not is_a_leaf) {
    int N_parts = gp_is.size();
    for (int i = 0; i < N_parts; i++) {
      int gp_i = gp_is[i];
      int gpp_i = gpp_is[i];
      scalar_type pp = p_part[i];
      if (not(e < last_leaf)) {
        int f = daughter[e];
        int g = son[e];
        // S event
        uq_sum += PS[e] * uq[gp_i][f] * uq[gpp_i][g] * pp + EPSILON;
        uq_sum += PS[e] * uq[gp_i][g] * uq[gpp_i][f] * pp + EPSILON;
      }
      // D event
      uq_sum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2) * pp + EPSILON;
      // T event
      for (int f = 0; f < last_branch; f++)
        if (not ancestral[e][f]) {
          uq_sum += uq[gp_i][e] * (PT[f] / (scalar_type) last_branch) * uq[gpp_i][f] * pp + EPSILON;
          uq_sum += uq[gpp_i][e] * (PT[f] / (scalar_type) last_branch) * uq[gp_i][f] * pp + EPSILON;
        }
    }
  }

  if (not(e < last_leaf)) {
    int f = daughter[e];
    int g = son[e];
    // SL event
    uq_sum += PS[e] * uq[i][f] * uE[g] + EPSILON;
    uq_sum += PS[e] * uq[i][g] * uE[f] + EPSILON;
  }

  // DL event
  uq_sum += PD[e] * (uq[i][e] * uE[e] * 2) + EPSILON;
  // TL event

  for (int f = 0; f < last_branch; f++)
    if (not ancestral[e][f]) {
      uq_sum += (PT[f] / (scalar_type) last_branch) * uq[i][f] * uE[e] + EPSILON;
      uq_sum += (PT[f] / (scalar_type) last_branch) * uE[f] * uq[i][e] + EPSILON;
    }

  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################

  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################

  std::stringstream estring;
  if (not(e < last_leaf)) estring << e; else estring << extant_species[e];
  std::string estr = estring.str();

  scalar_type uq_resum = 0;
  // S leaf and G leaf
  if (e < last_leaf and is_a_leaf and extant_species[e] == gid_sps[g_id]) {
    // present
    uq_resum += PS[e] * 1 + EPSILON;
    if (r * uq_sum < uq_resum) {
      register_leafu(e, last_event);
      return ale_pointer->set2name(ale_pointer->id_sets[g_id]) + branch_string + ":" + branch_length;
    }
  }
  // G internal
  if (not is_a_leaf) {
    int N_parts = gp_is.size();
    for (int i = 0; i < N_parts; i++) {
      int gp_i = gp_is[i];
      int gpp_i = gpp_is[i];
      scalar_type pp = p_part[i];
      if (not(e < last_leaf)) {
        int f = daughter[e];
        int g = son[e];
        // S event
        uq_resum += PS[e] * uq[gp_i][f] * uq[gpp_i][g] * pp + EPSILON;
        if (r * uq_sum < uq_resum) {
          register_Su(e, last_event);
          return "(" + sample_undated(f, gp_i, "S") + "," + sample_undated(g, gpp_i, "S") + ")." + estr +
                 branch_string + ":" + branch_length;
        }
        uq_resum += PS[e] * uq[gp_i][g] * uq[gpp_i][f] * pp + EPSILON;
        if (r * uq_sum < uq_resum) {
          register_Su(e, last_event);
          return "(" + sample_undated(g, gp_i, "S") + "," + sample_undated(f, gpp_i, "S") + ")." + estr +
                 branch_string + ":" + branch_length;
        }
      }
      // D event
      uq_resum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2) * pp + EPSILON;
      if (r * uq_sum < uq_resum) {
        register_D(e);
        return "(" + sample_undated(e, gp_i, "D") + "," + sample_undated(e, gpp_i, "D") + ").D@" + estr +
               branch_string + ":" + branch_length;
      }

      // T event
      for (int f = 0; f < last_branch; f++)
        if (not ancestral[e][f]) {
          std::stringstream fstring;
          if (not(f < last_leaf)) fstring << f; else fstring << extant_species[f];
          std::string fstr = fstring.str();

          uq_resum += uq[gp_i][e] * (PT[f] / (scalar_type) last_branch) * uq[gpp_i][f] * pp + EPSILON;
          if (r * uq_sum < uq_resum) {
            register_Tfrom(e);
            register_Tto(f);
            register_T_to_from(e, f);
            std::stringstream Ttoken;
            Ttoken << estr << ">" << fstr << "|" << ale_pointer->set2name(ale_pointer->id_sets[g_ids[gpp_i]]);
            Ttokens.push_back(Ttoken.str());

            return "(" + sample_undated(e, gp_i, "S") + "," + sample_undated(f, gpp_i, "T") + ").T@" + estr + "->" +
                   fstr + branch_string + ":" + branch_length;
          }
          uq_resum += uq[gpp_i][e] * (PT[f] / (scalar_type) last_branch) * uq[gp_i][f] * pp + EPSILON;
          if (r * uq_sum < uq_resum) {
            register_Tfrom(e);
            register_Tto(f);
            register_T_to_from(e, f);
            std::stringstream Ttoken;
            Ttoken << estr << ">" << fstr << "|" << ale_pointer->set2name(ale_pointer->id_sets[g_ids[gp_i]]);
            Ttokens.push_back(Ttoken.str());
            return "(" + sample_undated(e, gpp_i, "S") + "," + sample_undated(f, gp_i, "T") + ").T@" + estr + "->" +
                   fstr + branch_string + ":" + branch_length;
          }

        }
    }
  }
  if (not(e < last_leaf)) {
    int f = daughter[e];
    int g = son[e];
    // SL event
    uq_resum += PS[e] * uq[i][f] * uE[g] + EPSILON;
    if (r * uq_sum < uq_resum) {
      register_Su(e, last_event);
      register_L(g);
      return sample_undated(f, i, "S", "." + estr + branch_string);
    }
    uq_resum += PS[e] * uq[i][g] * uE[f] + EPSILON;
    if (r * uq_sum < uq_resum) {
      register_Su(e, last_event);
      register_L(f);
      return sample_undated(g, i, "S", "." + estr + branch_string);
    }
  }
  // DL event
  uq_resum += PD[e] * (uq[i][e] * uE[e] * 2) + EPSILON;
  if (r * uq_sum < uq_resum) {
    return sample_undated(e, i, "S", branch_string);
  }
  // TL event
  for (int f = 0; f < last_branch; f++)
    if (not ancestral[e][f]) {
      std::stringstream fstring;
      if (not(f < last_leaf)) fstring << f; else fstring << extant_species[f];
      std::string fstr = fstring.str();

      uq_resum += (PT[f] / (scalar_type) last_branch) * uq[i][f] * uE[e] + EPSILON;
      if (r * uq_sum < uq_resum) {
        register_Tfrom(e);
        register_Tto(f);
        register_T_to_from(e, f);

        //stringstream Ttoken;
        //Ttoken<<estr<<">"<<fstr<<"|"<<ale_pointer->set2name(ale_pointer->id_sets[g_id]);
        //Ttokens.push_back(Ttoken.str());

        register_L(e);
        return sample_undated(f, i, "T", ".T@" + estr + "->" + fstr + branch_string);
      }
      uq_resum += (PT[f] / (scalar_type) last_branch) * uE[f] * uq[i][e] + EPSILON;
      if (r * uq_sum < uq_resum) {
        return sample_undated(e, i, "S");
      }
    }
  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################
  std::cout << "sum error!" << std::endl;
  return "-!=-";
}

std::string exODT_model::counts_string_undated(scalar_type samples) {

  std::stringstream out;
  for (int e = 0; e < last_branch; e++) {
    bool isleaf = e < last_leaf;
    std::stringstream named_branch;
    if (e < last_leaf) named_branch << extant_species[e]; else named_branch << e;

    if (not isleaf)
      out << "S_internal_branch\t" << named_branch.str() << "\t"
          << branch_counts["Ds"][e] / samples << "\t"
          << branch_counts["Ts"][e] / samples << "\t"
          << branch_counts["Ls"][e] / samples << "\t"
          << branch_counts["Os"][e] / samples << "\t"
          //<< branch_counts["singleton"][e]/samples << "\t"
          << branch_counts["copies"][e] / samples << "\n";
    else
      out << "S_terminal_branch\t" << named_branch.str() << "\t"
          << branch_counts["Ds"][e] / samples << "\t"
          << branch_counts["Ts"][e] / samples << "\t"
          << branch_counts["Ls"][e] / samples << "\t"
          << branch_counts["Os"][e] / samples << "\t"
          //<< branch_counts["singleton"][e] << "\t"
          << branch_counts["copies"][e] / samples << "\n";

  }
  return out.str();
}

void exODT_model::register_Su(int e, std::string last_event) {
  MLRec_events["S"] += 1;
  if (e > -1) {
    int f = daughter[e];
    int g = son[e];
    if (last_event == "S" or last_event == "O") branch_counts["singleton"].at(e) += 1;
    branch_counts["copies"].at(e) += 1;
    branch_counts["count"].at(f) += 1;
    branch_counts["count"].at(g) += 1;
  }
}

void exODT_model::register_leafu(int e, std::string last_event) {
  if (e > -1) {
    branch_counts["copies"].at(e) += 1;
    if (last_event == "S" or last_event == "O") branch_counts["singleton"].at(e) += 1;
  }
  //MLRec_events["genes"]+=1;
}

void exODT_model::register_T_to_from(int e, int f) {
  T_to_from[e][f] += 1;
}

// FROM SAMPLE_SCALED.CPP
//
//consider reimplemtation for clarity!
//
//The general structure of the calculation, and lot of the code, is the same as p(ale) cf. model.cpp.
//(this could be made more clear)
string exODT_model::sample(bool max_rec)
{
  if (max_rec)
  {
    MLRec_events.clear();
    Ttokens.clear();
  }
  //scalar_type beta=1;
  scalar_type root_resum=0;
  for (int rank=0;rank<last_rank;rank++)
  {
    int n=time_slices[rank].size();
    for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
    {
      //scalar_type t=time_slice_times[rank][t_i];

      for (int branch_i=0;branch_i<n;branch_i++)
      {
        int e = time_slices[rank][branch_i];
        root_resum+=qvec[0][rank][t_i][e];
      }
      root_resum+=qvec[0][rank][t_i][alpha];
    }
  }
  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  scalar_type root_reresum=0;
  int root_e=-1;
  int root_rank=-1;
  int root_t_i=-1;

  int max_root_e=-1;
  int max_root_rank=-1;
  int max_root_t_i=-1;

  scalar_type max_resum=0;
  for (int rank=0;rank<last_rank;rank++)
  {
    int n=time_slices[rank].size();
    for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
    {
      //scalar_type t=time_slice_times[rank][t_i];

      for (int branch_i=0;branch_i<n;branch_i++)
      {
        int e = time_slices[rank][branch_i];
        root_reresum+=qvec[0][rank][t_i][e];
        if (max_resum<qvec[0][rank][t_i][e])
        {
          max_resum=qvec[0][rank][t_i][e];
          max_root_e=e;
          max_root_rank=rank;
          max_root_t_i=t_i;
        }
        if (r*root_resum<root_reresum and root_t_i==-1)
        {
          root_e=e;
          root_rank=rank;
          root_t_i=t_i;
        }
      }
      root_reresum+=qvec[0][rank][t_i][alpha];
      if (max_resum<qvec[0][rank][t_i][alpha])
      {
        max_resum=qvec[0][rank][t_i][alpha];
        max_root_e=alpha;
        max_root_rank=rank;
        max_root_t_i=t_i;
      }
      if (r*root_resum<root_reresum and root_t_i==-1)
      {
        root_e=alpha;
        root_rank=rank;
        root_t_i=t_i;
      }
    }
  }
  if (max_rec)
  {
    root_e=max_root_e;
    root_rank=max_root_rank;
    root_t_i=max_root_t_i;
  }
  register_O(root_e);
  /*
  if (root_e==-1)
    cout <<"c+O "<< time_slice_times[root_rank][root_t_i] << " " << 0 << endl;
  else
    cout <<"c+O "<< time_slice_times[root_rank][root_t_i] << " " << 1 << endl;
  */
  return sample(false,-1,root_t_i,root_rank,root_e,0,"","",max_rec)+";";
  //del-locs
}

//
//consider reimplemtation for clarity!
//
//The general structure of the calculation, and lot of the code, is the same as p(ale) cf. model.cpp.
//(this could be made more clear)
string exODT_model::sample(bool S_node,long int g_id,int t_i,scalar_type rank,int e,scalar_type branch_length,string branch_events, string transfer_token,bool max_rec)
{
  //cout << "iti " << t_i <<" " << rank<<" " << S_node << endl;
  // it could be nice to implemant a sampling temperature ?
  //scalar_type beta=1;
  stringstream topptmp;
  if (e==alpha)
    topptmp<<-1;
  else if (id_ranks[e]==0)
    topptmp<<extant_species[e];
  else
    topptmp<<id_ranks[e];

  auto ale = ale_pointer;
  bool is_a_leaf=false;

  int size = 0;
  boost::dynamic_bitset<> temp;
  //  if (g_id!=-1) { //We are not at the root bipartition
  temp = ale->id_sets.at( g_id );
  for (int i = 0; i < ale->Gamma_size + 1; ++i) {
    if ( temp[ i ] ) {
      size++;
    }
  }


  // if ((int)(ale->id_sets[g_id].size())==1)
  if (size == 1)
    is_a_leaf=true;
  //  }
  vector <long int> gp_ids;//del-loc
  vector <long int> gpp_ids;//del-loc
  //p_part is filled up  CCPs
  vector <scalar_type> p_part;//del-loc
  if (g_id!=-1)
    for (unordered_map< pair<long int, long int>,scalar_type> :: iterator kt = ale->Dip_counts[g_id].begin(); kt != ale->Dip_counts[g_id].end(); kt++)
    {
      pair<long int, long int> parts = (*kt).first;
      long int gp_id=parts.first;
      long int gpp_id=parts.second;
      gp_ids.push_back(gp_id);
      gpp_ids.push_back(gpp_id);
      if (ale->Bip_counts[g_id]<=scalar_parameter["min_bip_count"])
        p_part.push_back(0);
      else
      {
        p_part.push_back(ale->p_dip(g_id,gp_id,gpp_id));
      }
    }
  else
  {
    //root bipartition needs to be handled seperatly
    std::map<std::set<long int>,int> bip_parts;
    for (auto it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
    {
      long int gp_id=(*it).first;
      boost::dynamic_bitset<> gamma = ale->id_sets[gp_id];

      boost::dynamic_bitset<> not_gamma = ~gamma;
      not_gamma[0] = 0;

      /* for (auto i = 0; i < ale->nbint; ++i) {
   not_gamma[i] = 0;
   }
      /*
  for (set<int>::iterator st=ale->Gamma.begin();st!=ale->Gamma.end();st++)
  if (gamma.count(*st)==0)
  not_gamma.insert(*st);*/
      long int gpp_id = ale->set_ids[not_gamma];
      set <long int> parts;
      parts.insert(gp_id);
      parts.insert(gpp_id);
      bip_parts[parts]=1;
      //  gamma.clear();
      //  not_gamma.clear();
    }
    for (map<set<long int>,int> :: iterator kt = bip_parts.begin();kt!=bip_parts.end();kt++)
    {
      vector <long int> parts;
      for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) parts.push_back((*sit));
      long int gp_id=parts[0];
      long int gpp_id=parts[1];
      gp_ids.push_back(gp_id);
      gpp_ids.push_back(gpp_id);
      //if (ale->Bip_counts[gp_id]<=scalar_parameter["min_bip_count"] and not ale->Gamma_size<4)
      //  p_part.push_back(0);
      //else
      p_part.push_back(ale->p_bip(gp_id));
    }
    bip_parts.clear();
  }
  int N_parts=gp_ids.size();
  int n=time_slices[rank].size();
  //######################################################################################################################
  //######################################### INNER LOOP #################################################################
  //######################################################################################################################
  vector <step> sample_steps;
  vector <scalar_type> sample_ps;
  scalar_type resum=0;
  scalar_type t;

  //scalar_type t_to=time_slice_times[rank][t_i];

  //int rank_to=rank;
  //int t_i_to=t_i;
  bool set_S_node=false;
  // proceed a single "D" subslice
  if(t_i>0)
  {
    rank=rank;
    t_i-=1;
  }
    // at boundaries
  else if (rank>0)
  {
    if (S_node)
    {
      ;
    }
      //if e defines the time slice we have to look at speciaitons
    else if (e==time_slices[rank][n-1])
    {
      set_S_node=true;
    }
    else
    {
      rank-=1;
      t_i=time_slice_times[rank].size()-1;
    }
  }
  else
  {
    rank=-1;
    t_i=-1;
    if (is_a_leaf && extant_species[e]==gid_sps[g_id])
    {
      resum=1;
      sample_ps.push_back(1);
      step step;
      step.e=e;
      step.ep=-1;
      step.epp=-1;
      step.t=0;
      step.rank=0;
      step.g_id=g_id;
      step.gp_id=-1;
      step.gpp_id=-1;
      step.event="0";
      sample_steps.push_back(step);
    }//qvec[g_id+1][rank][t_i][e]=1;

  }
  if (rank>-1)
  {

    t=time_slice_times[rank][t_i];
    scalar_type tpdt;
    if ( t_i < (int)time_slice_times[rank].size()-1 )
      tpdt=time_slice_times[rank][t_i+1];
    else if (rank<last_rank-1)
      tpdt=time_slice_times[rank+1][0];
    else
      //top of root stem
      tpdt=t_begin[time_slices[rank][0]];

    //root
    scalar_type Delta_t=tpdt-t;
    //Delat_bar corresponds to \hat \sigma
    scalar_type ni=time_slices[rank].size();
    scalar_type delta_avg=scalar_parameter["delta_avg"];
    scalar_type tau_avg=scalar_parameter["tau_avg"];
    scalar_type lambda_avg=scalar_parameter["lambda_avg"];
    scalar_type sigma_hat=scalar_parameter["sigma_hat"];
    scalar_type H_hat=Ee[-1][t];


    if(e==alpha)
    {
      //boundaries for branch alpha virtual branch
      //boundary at present
      if (t==0)
      {
        resum+=0;
        sample_ps.push_back(0);
        step step;
        step.e=alpha;
        step.ep=-1;
        step.epp=-1;
        step.t=t;
        step.rank=rank;
        step.g_id=g_id;
        step.gp_id=-1;
        step.gpp_id=-1;
        step.event="0";
        sample_steps.push_back(step);
      }//qvec[g_id+1][rank][t_i][alpha]=0;
      //boundary between slice rank and rank-1 slice is trivial
      //trivial
      if (S_node )//and 0!?
      {
        resum=1;
        if(1)
        {
          sample_ps.push_back(1);
          step step;
          step.e=alpha;
          step.ep=-1;
          step.epp=-1;
          step.t=t;
          step.rank=rank;
          step.g_id=g_id;
          step.gp_id=-1;
          step.gpp_id=-1;
          step.event="0";
          sample_steps.push_back(step);
        }
        ;//qvec[g_id+1][rank][t_i][e]=qvec[g_id+1][rank][t_i][e];
      }
      else
      {
        //non-leaf directed partition
        if (not is_a_leaf)
          for (int i=0;i<N_parts;i++)
          {
            long int gp_id=gp_ids[i];
            long int gpp_id=gpp_ids[i];
            scalar_type pp=p_part[i];
            scalar_type Sb=2*sigma_hat*Delta_t*(qvec[gp_id+1][rank][t_i][alpha]*qvec[gpp_id+1][rank][t_i][alpha])*pp;
            resum+=Sb;
            if (1)
            {
              sample_ps.push_back(Sb);
              step step;
              step.e=-1;
              step.ep=alpha;
              step.epp=alpha;
              step.t=t;
              step.rank=rank;
              step.g_id=-1;
              step.gp_id=gp_id;
              step.gpp_id=gpp_id;
              step.event="Sb";
              sample_steps.push_back(step);
            }
          }


        for (int branch_i=0;branch_i<n;branch_i++)
        {
          int e = time_slices[rank][branch_i];
          scalar_type tau_e=vector_parameter["tau"][e];
          scalar_type TLb=tau_e*Delta_t*qvec[g_id+1][rank][t_i][e];
          //TL_bar EVENT
          resum+=TLb;
          if (1)
          {
            sample_ps.push_back(TLb);
            step step;
            step.e=e;
            step.ep=-1;
            step.epp=-1;
            step.t=t;
            step.rank=rank;
            step.g_id=g_id;
            step.gp_id=-1;
            step.gpp_id=-1;
            step.event="TLb";
            sample_steps.push_back(step);
          }
        }


        //0 EVENT
        scalar_type empty=(1+(delta_avg+tau_avg-lambda_avg-ni*sigma_hat-2*sigma_hat*H_hat)*Delta_t)*qvec[g_id+1][rank][t_i][alpha];
        resum+=empty;
        if (1)
        {
          sample_ps.push_back(empty);
          step step;
          step.e=alpha;
          step.ep=-1;
          step.epp=-1;
          step.t=t;
          step.rank=rank;
          step.g_id=g_id;
          step.gp_id=-1;
          step.gpp_id=-1;
          step.event="0";
          sample_steps.push_back(step);
        }
      }
    }
    else
    {


      //int e = time_slices[rank][branch_i];
      scalar_type Eet=Ee[e][t];
      scalar_type delta_e=vector_parameter["delta"][e];
      scalar_type lambda_e=vector_parameter["lambda"][e];

      if (S_node)
      {
        //boundaries for branch e
        //boundary at present
        if (t==0)
        {
          ;
        }
          //boundary between slice rank and rank-1
        else if (t_i==0 )
        {
          //terminating branch is last in time_slices and defines a represented speciation
          if (e==time_slices[rank][n-1] && rank>0)
          {
            int f=daughters[e][0];
            int g=daughters[e][1];
            scalar_type Eft=Ee[f][t];
            scalar_type Egt=Ee[g][t];

            //scalar_type q_sum=0;
            //qvec[g_id+1][rank][t_i][e]=0;

            scalar_type SL_fLg=qvec[g_id+1][rank][t_i][f]*Egt;
            scalar_type SL_Lfg=qvec[g_id+1][rank][t_i][g]*Eft;
            //SL EVENT
            resum+=SL_fLg;
            if(1)
            {
              sample_ps.push_back(SL_fLg);
              step step;
              step.e=f;
              step.ep=-1;
              step.epp=-1;
              step.t=t;
              step.rank=rank;
              step.g_id=g_id;
              step.gp_id=-1;
              step.gpp_id=-1;
              step.event="SL";
              sample_steps.push_back(step);
            }

            resum+=SL_Lfg;
            if(1)
            {
              sample_ps.push_back(SL_Lfg);
              step step;
              step.e=g;
              step.ep=-1;
              step.epp=-1;
              step.t=t;
              step.rank=rank;
              step.g_id=g_id;
              step.gp_id=-1;
              step.gpp_id=-1;
              step.event="SL";
              sample_steps.push_back(step);
            }


            //non-leaf directed partition
            if (not is_a_leaf)
              for (int i=0;i<N_parts;i++)
              {

                long int gp_id=gp_ids[i];
                long int gpp_id=gpp_ids[i];
                scalar_type pp=p_part[i];
                scalar_type S_pf_ppg=qvec[gp_id+1][rank][t_i][f]*qvec[gpp_id+1][rank][t_i][g]*pp;
                scalar_type S_ppf_pg=qvec[gpp_id+1][rank][t_i][f]*qvec[gp_id+1][rank][t_i][g]*pp;
                resum+=S_pf_ppg;
                if(1)
                {
                  sample_ps.push_back(S_pf_ppg);
                  step step;
                  step.e=-1;
                  step.ep=f;
                  step.epp=g;
                  step.t=t;
                  step.rank=rank;
                  step.g_id=-1;
                  step.gp_id=gp_id;
                  step.gpp_id=gpp_id;
                  step.event="S";
                  sample_steps.push_back(step);
                }

                resum+=S_ppf_pg;
                if(1)
                {
                  sample_ps.push_back(S_ppf_pg);
                  step step;
                  step.e=-1;
                  step.ep=g;
                  step.epp=f;
                  step.t=t;
                  step.rank=rank;
                  step.g_id=-1;
                  step.gp_id=gp_id;
                  step.gpp_id=gpp_id;
                  step.event="S";
                  sample_steps.push_back(step);
                }

              }
          }
            //branches that cross to next time slice
          else
          {

            //trivial
            resum=1;
            if(1)
            {
              sample_ps.push_back(1);
              step step;
              step.e=e;
              step.ep=-1;
              step.epp=-1;
              step.t=t;
              step.rank=rank;
              step.g_id=g_id;
              step.gp_id=-1;
              step.gpp_id=-1;
              step.event="0";
              sample_steps.push_back(step);
            }
            ;//qvec[g_id+1][rank][t_i][e]=qvec[g_id+1][rank][t_i][e];
          }
        }
      }
        //boundaries for branch e.
      else
      {



        //non-leaf directed partition
        if (not is_a_leaf)
          for (int i=0;i<N_parts;i++)
          {


            long int gp_id=gp_ids[i];
            long int gpp_id=gpp_ids[i];
            scalar_type pp=p_part[i];
            scalar_type qpe=qvec[gp_id+1][rank][t_i][e];
            scalar_type qppe=qvec[gpp_id+1][rank][t_i][e];
            scalar_type Sb_pa_ppe= sigma_hat*Delta_t*qvec[gp_id+1][rank][t_i][alpha]*qppe*pp;
            scalar_type Sb_pe_ppa= sigma_hat*Delta_t*qpe*qvec[gpp_id+1][rank][t_i][alpha]*pp;

            //S_bar EVENT
            resum+=Sb_pa_ppe;
            if(1)
            {
              sample_ps.push_back(Sb_pa_ppe);
              step step;
              step.e=-1;
              step.ep=alpha;
              step.epp=e;
              step.t=t;
              step.rank=rank;
              step.g_id=-1;
              step.gp_id=gp_id;
              step.gpp_id=gpp_id;
              step.event="Sb";
              sample_steps.push_back(step);
            }


            resum+=Sb_pe_ppa;
            if(1)
            {
              sample_ps.push_back(Sb_pe_ppa);
              step step;
              step.e=-1;
              step.ep=e;
              step.epp=alpha;
              step.t=t;
              step.rank=rank;
              step.g_id=-1;
              step.gp_id=gp_id;
              step.gpp_id=gpp_id;
              step.event="Sb";

              sample_steps.push_back(step);
            }

            scalar_type D=2*delta_e*Delta_t*qpe*qppe*pp;
            resum+=D;
            if(1)
            {
              sample_ps.push_back(D);
              step step;
              step.e=-1;
              step.ep=e;
              step.epp=e;
              step.t=t;
              step.rank=rank;
              step.g_id=-1;
              step.gp_id=gp_id;
              step.gpp_id=gpp_id;
              step.event="D";

              sample_steps.push_back(step);
            }
            //D EVENT
            //q_sum+= D;
            //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=2*p_delta_e*qvec[gp_id+1][rank][t_i][e]*qvec[gpp_id+1][rank][t_i][e];
            //D.

          }
        scalar_type SLb=sigma_hat*Delta_t*Eet*qvec[g_id+1][rank][t_i][alpha];
        //SL_bar EVENT
        resum+=SLb;
        if(1)
        {
          sample_ps.push_back(SLb);
          step step;
          step.e=alpha;
          step.ep=-1;
          step.epp=-1;
          step.t=t;
          step.rank=rank;
          step.g_id=g_id;
          step.gp_id=-1;
          step.gpp_id=-1;
          step.event="SLb";
          sample_steps.push_back(step);
        }
        //q_suml+=SLb;

        //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=p_Delta_bgar*Eet*qvec[g_id+1][rank][t_i][alpha];
        //SL_bar.


        scalar_type empty=(1+(2*delta_e*Eet-sigma_hat*H_hat-delta_e-lambda_e)*Delta_t)*qvec[g_id+1][rank][t_i][e];
        //0 EVENT
        resum+=empty;
        if(1)
        {
          sample_ps.push_back(empty);
          step step;
          step.e=e;
          step.ep=-1;
          step.epp=-1;
          step.t=t;
          step.rank=rank;
          step.g_id=g_id;
          step.gp_id=-1;
          step.gpp_id=-1;
          step.event="0";
          sample_steps.push_back(step);
        }
        //q_sum+=empty;

        //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]=Get*qvec[g_id+1][rank][t_i][e];
        //0.

        //qvec[g_id+1][tpdt_rank][tpdt_t_i][e]+=q_sum;

        //events within slice rank at time t on branch e.
      }
    }
  }
  //######################################################################################################################
  //#########################################INNNER LOOP##################################################################
  //######################################################################################################################
  gp_ids.clear();
  gpp_ids.clear();
  p_part.clear();
  if (S_node)
  {
    rank-=1;
    t_i=time_slice_times[rank].size()-1;
    S_node=false;
  }
  if (set_S_node)
  {
    S_node=true;
  }
  int step_i=-1;
  scalar_type reresum=0;
  scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  scalar_type max_resum=0;
  int max_i=0;
  for (int i=0;i<(int)sample_ps.size();i++)
  {
    //cout <<i << " " << sample_ps[i] << " " << sample_steps[i].event << endl;;
    if (max_resum<sample_ps[i])
    {
      max_resum=sample_ps[i];
      max_i=i;
    }
    reresum+=sample_ps[i];
    if (r*resum<reresum )
    {
      step_i=i;
      if (not max_rec)
        break;
    }
  }
  if (max_rec)
    step_i=max_i;
  step back_step=sample_steps.at(step_i);
  sample_steps.clear();
  sample_ps.clear();
  /*
    if (back_step.e==-1)
      cout  <<  back_step.event << "\t" << back_step.rank << "\tid_rank:" << id_ranks[back_step.ep] << "\text_sp:" << extant_species[back_step.ep] << "\te:" << back_step.ep<< "\tg_id:" << ale_pointer->set2name(ale_pointer->id_sets[g_id]) << "\t" << back_step.t << " t_i:" << t_i << " " << time_slice_times[back_step.rank][t_i]<< endl;
    else
      cout  <<  back_step.event << "\t" << back_step.rank << "\tid_rank:" << id_ranks[back_step.e] << "\text_sp:" << extant_species[back_step.e] << "\te:" << back_step.e<< "\tg_id:" << ale_pointer->set2name(ale_pointer->id_sets[g_id]) << "\t" << back_step.t << " t_i:" << t_i  << " " << time_slice_times[back_step.rank][t_i]<<  endl;
  */
  stringstream toptmp;
  if (back_step.e==alpha)
    toptmp<<-1;
  else if (id_ranks[back_step.e]==0)
    toptmp<<extant_species[back_step.e];
  else
    toptmp<<id_ranks[back_step.e];
  //cout << branch_length << " +  " << t << " + " << back_step.t << endl;
  scalar_type new_branch_length=-1;//branch_length+t-back_step.t;


  size = 0;
  //  if (g_id!=-1) { //We are not at the root bipartition

  temp = ale_pointer->id_sets[g_id];
  for (int i = 0; i < ale_pointer->Gamma_size + 1; ++i) {
    if ( temp[i] ) {
      size++;
    }
  }
  //   }


  if (ale_pointer->Bip_counts.count(g_id) and ale_pointer->Bip_counts[g_id]>0)
  {
    new_branch_length=max(ale_pointer->Bip_bls[g_id]/ale_pointer->Bip_counts[g_id],(scalar_type)scalar_parameter["min_branch_lenghts"]);
  }
  else
  {
    new_branch_length=max(ale_pointer->Bip_bls[g_id]/ale_pointer->observations,(scalar_type)scalar_parameter["min_branch_lenghts"]);

  }

  if (back_step.t==0 and size == 1 and e!=-1)
  {
    //cout <<"c+P "<< back_step.t << " " << 0 << endl;

    register_leaf(e);
    stringstream branch_string;
    if (scalar_parameter["leaf_events"]==1) branch_string<<branch_events;
    branch_string <<":"<<new_branch_length;

    gid_events[g_id].push_back(">PRESENT");
    gid_times[g_id].push_back(0);
    gid_branches[g_id].push_back(e);
    gid_gidp[g_id].push_back(g_id);
    gid_gidpp[g_id].push_back(g_id);

    return ale_pointer->set2name(ale_pointer->id_sets[g_id])+branch_string.str();
  }

  if (back_step.event=="D" or back_step.event=="Tb" or back_step.event=="S" or back_step.event=="Sb")
  {

    stringstream transfer_token_stream;
    transfer_token_stream<<"";
    stringstream branch_string;
    if (back_step.event=="S")
    {
      //cout <<"c+S "<< back_step.t << " " << 1 << endl;

      register_S(e);

      gid_events[g_id].push_back(">S");
      gid_times[g_id].push_back(t);
      gid_branches[g_id].push_back(e);
      gid_gidp[g_id].push_back(back_step.gp_id);
      gid_gidpp[g_id].push_back(back_step.gpp_id);

      gid_events[back_step.gp_id].push_back("S<");
      gid_times[back_step.gp_id].push_back(t);
      gid_branches[back_step.gp_id].push_back(back_step.ep);
      gid_gidp[back_step.gp_id].push_back(back_step.gp_id);
      gid_gidpp[back_step.gp_id].push_back(back_step.gp_id);

      gid_events[back_step.gpp_id].push_back("S<");
      gid_times[back_step.gpp_id].push_back(t);
      gid_branches[back_step.gpp_id].push_back(back_step.epp);
      gid_gidp[back_step.gpp_id].push_back(back_step.gpp_id);
      gid_gidpp[back_step.gpp_id].push_back(back_step.gpp_id);

      branch_string<< branch_events
                   <<"."<<id_ranks[e]<<":"<<max(new_branch_length,(scalar_type)0.0);
    }
    else
    {
      if (back_step.event=="Tb")
      {
        //cout <<"c+Tb "<< back_step.t << " " << 1 << endl;

        int this_e,this_gid;

        gid_events[g_id].push_back(">T");
        gid_times[g_id].push_back(t);
        gid_branches[g_id].push_back(e);


        if (back_step.ep==alpha)
        {
          this_e=back_step.epp;
          this_gid=back_step.gpp_id;

          gid_gidp[g_id].push_back(back_step.gpp_id);
          gid_gidpp[g_id].push_back(back_step.gp_id);

          gid_events[back_step.gp_id].push_back("Tto<");
          gid_times[back_step.gp_id].push_back(t);
          gid_branches[back_step.gp_id].push_back(back_step.ep);
          gid_gidp[back_step.gp_id].push_back(back_step.gp_id);
          gid_gidpp[back_step.gp_id].push_back(back_step.gp_id);

          gid_events[back_step.gpp_id].push_back("Tfrom<");
          gid_times[back_step.gpp_id].push_back(t);
          gid_branches[back_step.gpp_id].push_back(back_step.epp);
          gid_gidp[back_step.gpp_id].push_back(back_step.gpp_id);
          gid_gidpp[back_step.gpp_id].push_back(back_step.gpp_id);

        }
        else
        {
          this_e=back_step.ep;
          this_gid=back_step.gp_id;

          gid_gidp[g_id].push_back(back_step.gp_id);
          gid_gidpp[g_id].push_back(back_step.gpp_id);

          gid_events[back_step.gp_id].push_back("Tfrom<");
          gid_times[back_step.gp_id].push_back(t);
          gid_branches[back_step.gp_id].push_back(back_step.ep);
          gid_gidp[back_step.gp_id].push_back(back_step.gp_id);
          gid_gidpp[back_step.gp_id].push_back(back_step.gp_id);

          gid_events[back_step.gpp_id].push_back("Tto<");
          gid_times[back_step.gpp_id].push_back(t);
          gid_branches[back_step.gpp_id].push_back(back_step.epp);
          gid_gidp[back_step.gpp_id].push_back(back_step.gpp_id);
          gid_gidpp[back_step.gpp_id].push_back(back_step.gpp_id);

        }

        stringstream named_branch;
        if (this_e==alpha)
          named_branch<<-1;
        else if (id_ranks[this_e]==0)
          named_branch<<extant_species[this_e];
        else
          named_branch<<id_ranks[this_e];
        // Tto
        register_Tto(this_e);


        stringstream tmp;
        tmp<<back_step.rank<<"|"<<t<<"|"<<named_branch.str()<<"|"<<this_gid;
        register_Ttoken(transfer_token+">"+tmp.str());
        // Tto

        branch_string<< branch_events<<back_step.event<<"@"<<back_step.rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0);
      }
      else if (back_step.event=="Sb")
      {
        int this_e;
        if (back_step.ep==alpha)
        {
          this_e=back_step.epp;

          gid_gidp[g_id].push_back(back_step.gpp_id);
          gid_gidpp[g_id].push_back(back_step.gp_id);

          gid_events[back_step.gp_id].push_back("Sto<");
          gid_times[back_step.gp_id].push_back(t);
          gid_branches[back_step.gp_id].push_back(back_step.ep);
          gid_gidp[back_step.gp_id].push_back(back_step.gp_id);
          gid_gidpp[back_step.gp_id].push_back(back_step.gp_id);

          gid_events[back_step.gpp_id].push_back("Sfrom<");
          gid_times[back_step.gpp_id].push_back(t);
          gid_branches[back_step.gpp_id].push_back(back_step.epp);
          gid_gidp[back_step.gpp_id].push_back(back_step.gpp_id);
          gid_gidpp[back_step.gpp_id].push_back(back_step.gpp_id);

        }
        else
        {
          this_e=back_step.ep;

          gid_gidp[g_id].push_back(back_step.gp_id);
          gid_gidpp[g_id].push_back(back_step.gpp_id);

          gid_events[back_step.gp_id].push_back("Sfrom<");
          gid_times[back_step.gp_id].push_back(t);
          gid_branches[back_step.gp_id].push_back(back_step.ep);
          gid_gidp[back_step.gp_id].push_back(back_step.gp_id);
          gid_gidpp[back_step.gp_id].push_back(back_step.gp_id);

          gid_events[back_step.gpp_id].push_back("Sto<");
          gid_times[back_step.gpp_id].push_back(t);
          gid_branches[back_step.gpp_id].push_back(back_step.epp);
          gid_gidp[back_step.gpp_id].push_back(back_step.gpp_id);
          gid_gidpp[back_step.gpp_id].push_back(back_step.gpp_id);

        }
        // Tfrom
        register_Tfrom(this_e);

        gid_events[g_id].push_back(">Sfrom");
        gid_times[g_id].push_back(t);
        gid_branches[g_id].push_back(this_e);

        // Tfrom
        stringstream named_branch;
        if (this_e==alpha)
          named_branch<<-1;
        else if (id_ranks[this_e]==0)
          named_branch<<extant_species[this_e];
        else
          named_branch<<id_ranks[this_e];
        if  (transfer_token!="")
          transfer_token_stream<< transfer_token;
        else
          transfer_token_stream<<"T|"<< rank<<"|"<<t<<"|"<<named_branch.str()<<"|"<<g_id;

        branch_string<< branch_events<<"T@"<<rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0);
      }
      else
      {
        //cout <<"c+D "<< back_step.t << " " << 1 << endl;
        register_D(e);
        stringstream Dtoken_stream;
        stringstream named_branch;
        if (e==alpha)
          named_branch<<-1;
        else if (id_ranks[e]==0)
          named_branch<<extant_species[e];
        else
          named_branch<<id_ranks[e];

        gid_events[g_id].push_back(">D");
        gid_times[g_id].push_back(t);
        gid_branches[g_id].push_back(e);
        gid_gidp[g_id].push_back(back_step.gp_id);
        gid_gidpp[g_id].push_back(back_step.gpp_id);

        gid_events[back_step.gp_id].push_back("Dto<");
        gid_times[back_step.gp_id].push_back(t);
        gid_branches[back_step.gp_id].push_back(back_step.ep);
        gid_gidp[back_step.gp_id].push_back(back_step.gp_id);
        gid_gidpp[back_step.gp_id].push_back(back_step.gp_id);

        gid_events[back_step.gpp_id].push_back("Dfrom<");
        gid_times[back_step.gpp_id].push_back(t);
        gid_branches[back_step.gpp_id].push_back(back_step.epp);
        gid_gidp[back_step.gpp_id].push_back(back_step.gpp_id);
        gid_gidpp[back_step.gpp_id].push_back(back_step.gpp_id);


        Dtoken_stream    << "D|" << rank << "|" <<named_branch.str() << "|"<< g_id;
        register_Ttoken(Dtoken_stream.str());

        branch_string<< branch_events<<back_step.event<<"@"<<rank<<"|"<<named_branch.str()<<":"<<max(new_branch_length,(scalar_type)0.0);
      }
    }
    if ( transfer_token_stream.str()=="")
      return "("+
             sample(S_node,back_step.gp_id, t_i+(back_step.event=="S"
             ),rank,back_step.ep,0,"",transfer_token,max_rec)
             +","+
             sample(S_node,back_step.gpp_id, t_i+(back_step.event=="S"
             ),rank,back_step.epp,0,"",transfer_token,max_rec)
             +")"+branch_string.str();
    else
    if(back_step.ep==alpha)
      return "("+
             sample(S_node,back_step.gp_id, t_i+(back_step.event=="S"
             ),rank,back_step.ep,0,"",transfer_token_stream.str(),max_rec)
             +","+
             sample(S_node,back_step.gpp_id, t_i+(back_step.event=="S"
             ),rank,back_step.epp,0,"",transfer_token,max_rec)
             +")"+branch_string.str();
    else
      return "("+
             sample(S_node,back_step.gp_id, t_i+(back_step.event=="S"
             ),rank,back_step.ep,0,"",transfer_token)
             +","+
             sample(S_node,back_step.gpp_id, t_i+(back_step.event=="S"
             ),rank,back_step.epp,0,"",transfer_token_stream.str(),max_rec)
             +")"+branch_string.str();

  }
  else if ( back_step.event=="TLb" or back_step.event=="SL" or back_step.event=="SLb" or back_step.event=="0")
  {

    stringstream branch_string;
    stringstream transfer_token_stream;
    transfer_token_stream <<"";

    branch_string<< branch_events;
    if (back_step.event!="0")
    {
      if (back_step.event=="SL")
      {
        t_i=time_slice_times[rank].size()-1;
        register_S(e);
        int f=daughters[e][0];
        int g=daughters[e][1];
        if (back_step.e==f)
          register_L(g);
        else
          register_L(f);
        branch_string<<"."
                     <<id_ranks[e];

        gid_events[g_id].push_back("SL");
        gid_times[g_id].push_back(t);
        gid_branches[g_id].push_back(back_step.e);
        gid_gidp[g_id].push_back(g_id);
        gid_gidpp[g_id].push_back(g_id);


      }
      else
      {
        if (back_step.event=="TLb")
        {
          //cout <<"c+TLb "<< back_step.t << " " << (back_step.e!=-1) << endl;
          register_Tto(back_step.e);

          gid_events[g_id].push_back("TL");
          gid_times[g_id].push_back(t);
          gid_branches[g_id].push_back(back_step.e);
          gid_gidp[g_id].push_back(g_id);
          gid_gidpp[g_id].push_back(g_id);

          stringstream tmp;

          stringstream named_branch;
          if (back_step.e==alpha)
            named_branch<<-1;
          else if (id_ranks[back_step.e]==0)
            named_branch<<extant_species[back_step.e];
          else
            named_branch<<id_ranks[back_step.e];

          tmp<<back_step.rank<<"|"<<t<<"|"<< named_branch.str()<<"|"<<g_id;
          register_Ttoken(transfer_token+">"+tmp.str());
          transfer_token="";

          branch_string<<""
                       <<"@"<<back_step.rank<<"|"<<named_branch.str();
        }
        else  if (back_step.event=="SLb")
        {
          //cout <<"c+TLb "<< back_step.t << " " << -1 << endl;
          register_L(e);
          register_Tfrom(e);

          gid_events[g_id].push_back("SLb");
          gid_times[g_id].push_back(t);
          gid_branches[g_id].push_back(e);
          gid_gidp[g_id].push_back(g_id);
          gid_gidpp[g_id].push_back(g_id);

          stringstream named_branch;
          if (e==alpha)
            named_branch<<-1;
          else if (id_ranks[e]==0)
            named_branch<<extant_species[e];
          else
            named_branch<<id_ranks[e];

          transfer_token_stream<<"T|"<< rank<<"|"<<t<<"|"<<named_branch.str()<<"|"<<g_id;

          branch_string<<".T"
                       <<"@"<<rank<<"|"<<named_branch.str();
        }
      }
    }
    if (transfer_token_stream.str()=="")
    {
      return sample(S_node,back_step.g_id, t_i+(back_step.event=="SL"),rank, back_step.e, new_branch_length, branch_string.str(),transfer_token,max_rec);
    }
    else
    {
      return sample(S_node,back_step.g_id, t_i+(back_step.event=="SL"),rank, back_step.e, new_branch_length, branch_string.str(),transfer_token_stream.str(),max_rec);
    }
  }
  else
  {
    cout << "error "  <<endl;
    cout << " g_id " << g_id;
    cout << " t " << t;
    cout << " e " << e;
    cout << " l " << branch_length;
    cout << " str " << branch_events;
    cout << endl;
    cout << ale_pointer->constructor_string <<endl;
//    signal=-11;
  }
  return "error";
}

// FROM TRACEBACK_SCALED.CPP

//The current lowmem=true method uses sample(true) cf. sample.cpp.
//The general structure of the calculation, and lot of the code, is the same as p(ale) cf. model.cpp.

//used by sample() consider moving to sample.cpp
void exODT_model::register_O(int e)
{
  if (e>-1) branch_counts["count"].at(e)+=1;
  if (e>-1) branch_counts["Os"].at(e)+=1;
}
void exODT_model::register_D(int e)
{
  MLRec_events["D"]+=1;
  if (e>-1) branch_counts["Ds"].at(e)+=1;
}

void exODT_model::register_Tto(int e)
{
  MLRec_events["T"]+=1;
  if (e>-1) branch_counts["Ts"].at(e)+=1;
}

void exODT_model::register_Tfrom(int e)
{
  if (e>-1) branch_counts["Tfroms"].at(e)+=1;
}

void exODT_model::register_L(int e)
{
  MLRec_events["L"]+=1;
  if (e>-1) branch_counts["Ls"].at(e)+=1;
}
void exODT_model::register_S(int e)
{
  MLRec_events["S"]+=1;
  if (e>-1)
  {
    int f=daughters[e][0];
    int g=daughters[e][1];
    branch_counts["copies"].at(e)+=1;
    branch_counts["count"].at(f)+=1;
    branch_counts["count"].at(g)+=1;
  }
}

void exODT_model::register_leaf(int e)
{
  if (e>-1) branch_counts["copies"].at(e)+=1;
  //MLRec_events["genes"]+=1;
}

void exODT_model::register_Ttoken(string token)
{
  Ttokens.push_back(token);
}

//ad hoc function should be moved to a future exODT_util.cpp
string exODT_model::vertical_string(long int g_id, string ancestral_string,scalar_type t_0)
{
  stringstream event_stream;
  if (t_0==-1) t_0=gid_times[g_id][0];
  for (int i = 0; i < (int)gid_branches[g_id].size(); i++)
  {
    int branch = gid_branches[g_id][i];
    stringstream named_branch;
    if (branch==alpha)
      named_branch<<-1;
    else if (id_ranks[branch]==0)
    {
      named_branch<<extant_species[branch];
    }
    else
      named_branch<<id_ranks[branch];

    event_stream<<gid_events[g_id][i]
                <<"@"<<named_branch.str()
                //<<"*"<<gid_times[g_id][i]
                <<"|";
  }
  int last_i=(int)gid_branches[g_id].size()-1;
  if (not (gid_events[g_id][last_i]==">PRESENT" or gid_events[g_id][last_i]==">S" or ( gid_events[g_id][last_i]==">Sfrom" and gid_gidp[g_id][last_i]!=gid_gidpp[g_id][last_i]) or gid_events[g_id][last_i]==">D"))
  {
    return vertical_string(gid_gidp[g_id][last_i],event_stream.str(),t_0);
  }
  else
  {
    scalar_type bnorm;//=ale_pointer->Bip_counts[g_id];
    if (ale_pointer->Bip_counts.count(g_id)==0 )
      bnorm=ale_pointer->observations;
    else
      bnorm=ale_pointer->Bip_counts[g_id];
    if (bnorm==0)
      bnorm=ale_pointer->observations;

    event_stream<<" "<<t_0<< " " << gid_times[g_id][last_i] << " " << ale_pointer->Bip_bls[g_id]/bnorm <<"\t"<< gid_gidp[g_id][last_i] << "\t" << gid_gidpp[g_id][last_i];
    return ancestral_string+event_stream.str();
  }
}

