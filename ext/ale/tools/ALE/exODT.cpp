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
  // This, labels the bootstrap.
  // Each branch is edited with a property "ID".
  for (auto it = node_ids.begin(); it != node_ids.end(); it++) {
    auto node = (*it).first;
    int branch = (*it).second;
    std::stringstream out;
    std::stringstream out1;
    std::stringstream out2;
    out1 << 0.0;
    out2 << 0.0;
    int rank = branch;
    out << rank;
    if (S->hasFather(node) and S->getEdgeToFather(node))
      S->getEdgeToFather(node)->setProperty("ID", BppString(out.str()));
  }

  // Save the resulting speciestree with the "ID" label at each branch.
  string_parameter["S_with_ranks"] = IO::PhyloTreeToNewick(*S, false, "ID");

  // Initialize the ancestor matrix (nb_branches x bn_branches) with zeros.
  // Value take one if node has parent's relationship with an other.
  ancestral.clear();
  ancestral.resize(last_branch);
  for (int e = 0; e < last_branch; e++) {
    ancestral[e].resize(last_branch);
    ancestors.push_back(std::vector<int>());
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
    while (S->hasFather(node)) { //While we are not at the root.
      std::stringstream name_to;
      int f = node_ids[node];
      if (f < last_leaf)
        name_to << node_name[node];
      else
        name_to << f;

      node = S->getFather(node);
      if (not ancestral[edge][f])
        ancestors[edge].push_back(f);

      ancestral[edge][f] = 1;
    }
    std::stringstream name_to;
    int f = node_ids[node];
    name_to << f;
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
    if (i > 0) // There should be no need for this loop at the first iteration, because then it leaves mPTE_ancestral_correction at 0.
    {
      for (int edge = 0; edge < last_branch; edge++) {
        mPTE_ancestral_correction[edge] = 0;
        for (auto it = ancestors[edge].begin(); it != ancestors[edge].end(); it++) {
          int f = (*it);
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

  //directed partitions and their sizes
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
  
exODT_model::~exODT_model() {
    extant_species.clear();
    scalar_parameter.clear();
    for (std::map<std::string, std::vector<scalar_type> >::iterator it = vector_parameter.begin();
         it != vector_parameter.end(); it++)//del_loc
      (*it).second.clear();
    vector_parameter.clear();
    string_parameter.clear();
    node_ids.clear();
//      delete S;
    gid_sps.clear();

  }

