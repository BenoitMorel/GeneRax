// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 14/12/16.
//

// Include Bpp
#include <Bpp/BppString.h>

// Include Treerecs
#include "SpeciesGeneMapper.h"
#include <ale/tools/IO/IO.h>

bool SpeciesGeneMapper::case_sensitive_ = DEFAULT_CASE_SENSITIVE_MAPPING;

std::string
SpeciesGeneMapper::get_specie_from_gene_name_(const std::string &genename, const std::string& separator, const bool before) {
  auto genename_split = Utils::splitString(genename, separator.c_str());

  std::string species_name;

  if(before) {
    for(std::size_t i = 0; i < genename_split.size() - 1; ++i) {
      species_name += genename_split.at(i);
      if(i < genename_split.size() - 2) species_name += separator;
    }

  } else {
    for(std::size_t i = 1; i < genename_split.size(); ++i) {
      species_name += genename_split.at(i);
      if(i < (genename_split.size() - 1))
        species_name += separator;
    }
  }

  return species_name;
}

bool SpeciesGeneMapper::mappableWithNhxTags(
    const bpp::PhyloTree &genetree, const bpp::PhyloTree &speciestree
    , const std::string &tag
) {
  auto genetree_leaves = genetree.getAllLeaves();
  for(auto leaf: genetree_leaves) {
    if(not leaf->hasProperty(tag)) {
      return false;
    }
  }
  return true;
}

SpeciesGeneMap SpeciesGeneMapper::mapWithNhxTags(
    const bpp::PhyloTree &genetree, const bpp::PhyloTree &speciestree
    , const std::string &tag
) {
  auto gene_leaves = genetree.getAllLeaves();
  auto species_leaves = speciestree.getAllLeaves();
  SpeciesGeneMap map;
  std::list<Node> genes_without_tag;

  // Transform speciestree_leaves into a umap with the node name as key.
  std::unordered_map<std::string, Node> species_leaves_dictionary;
  for(auto snode: species_leaves){
    std::string sname = snode->getName();
    if(not case_sensitive_)
      std::transform(sname.begin(), sname.end(), sname.begin(), ::tolower);
    species_leaves_dictionary[sname] = snode;
  }

  for(auto gene_leaf: gene_leaves) {
    if(gene_leaf->hasProperty(tag)) {
      std::string species_name =
          dynamic_cast<bpp::BppString*>(gene_leaf->getProperty(tag))->toSTL();

      if(not case_sensitive_)
        std::transform(species_name.begin(), species_name.end(), species_name.begin(), ::tolower);

      if(species_leaves_dictionary.find(species_name) ==
         species_leaves_dictionary.end()) {
        std::cerr << "Error during species<>gene mapping: species "
                  << species_name;
        if(gene_leaf->hasName()) {
          std::cerr << " (gene " << gene_leaf->getName() << ")";
        }
        std::cerr << " does not exist in speciestree." << std::endl;
        exit(EXIT_FAILURE);
      }
      Node species_leaf = species_leaves_dictionary[species_name];
      map.addGene(species_leaf, gene_leaf);
    } else {
      // Save a list of genes with missing tags
      genes_without_tag.push_back(gene_leaf);
    }
  }

  if(not genes_without_tag.empty()) {
    std::cerr << "Error during species<>node mapping, somes genes ";
    std::cerr << "have no associated tag in NHX to get a species:" << std::endl;
    for(auto err_genes: genes_without_tag) {
      if(err_genes->hasName()) {
        std::cerr << "* " << err_genes->getName() << std::endl;
      } else {
        std::cerr << "* anonymous gene leaf" << std::endl;
      }
    }
  }

  return map;
}

GeneMap<std::string, std::string> SpeciesGeneMapper::mapWithNhxTags(
    const std::string &genetree, const std::string &speciestree
    , const std::string &tag
) {
  auto gene_tree = IO::nhxToPhyloTree(genetree);
  auto species_tree = IO::newickToPhyloTree(speciestree);

  auto map = SpeciesGeneMapper::mapWithNhxTags(*gene_tree, *species_tree, tag);

  return nodeMapsToStrings(map, *gene_tree, *species_tree, true);
}

SpeciesGeneMap
SpeciesGeneMapper::mapWithGenesNames(
    const bpp::PhyloTree &genetree,
    const bpp::PhyloTree &speciestree,
    std::string separator, const bool before) {
  SpeciesGeneMap ret;

  auto species_leaves = speciestree.getAllLeaves();

  for(std::shared_ptr<bpp::PhyloNode>& genenode: genetree.getAllLeaves()) {
    if(genenode->hasName()) {
      //If the gene node has a name we can looking for a species name inside.
      std::string gene_name = genenode->getName();
      std::string species_name = get_specie_from_gene_name_(gene_name, separator, before);
      std::shared_ptr<bpp::PhyloNode> species_node = PhyloTreeToolBox::getNodeFromName(species_leaves, species_name, case_sensitive_);
      if(species_node == nullptr){
        std::cerr << "There is no species associated with the gene, named: " << gene_name << std::endl << std::endl;
      } else {
        ret.addGene(species_node, genenode);
      }
    }
  }
  return ret;
}

GeneMap<std::string, std::string> SpeciesGeneMapper::mapWithGenesNames(
    const std::string &genetree,
    const std::string &speciestree,
    std::string separator,
    const bool before) {
  auto gene_tree = IO::newickToPhyloTree(genetree);
  auto species_tree = IO::newickToPhyloTree(speciestree);

  auto map = mapWithGenesNames(*gene_tree, *species_tree, separator, before);

  return nodeMapsToStrings(map, *gene_tree, *species_tree, true);
}

SpeciesGeneMap
SpeciesGeneMapper::mapWithTrees(const bpp::PhyloTree &genetree, const bpp::PhyloTree &speciestree) {

  // Initialization of the object returned by the function
  SpeciesGeneMap genemap_returned;

  // Get all nodes of each tree.
  std::vector<std::shared_ptr<bpp::PhyloNode>> species_nodes = speciestree.getAllLeaves();
  std::vector<std::shared_ptr<bpp::PhyloNode>> genes_nodes = genetree.getAllLeaves();
  std::list<Node> unmapped_gene_nodes;

  for( auto& gene_node: genes_nodes ){ // For each gene node
    if( gene_node->hasName() ) { // if the node has a name

      // Init a list of species candidates which can be associated with the gene name.
      std::list<std::pair<std::shared_ptr<bpp::PhyloNode>, bool>> species_candidates;
      // list of species candidates which has a match with a name in the gene name.
      //species_candidates.resize(species_nodes.size());

      // Get the name of the gene node.
      std::string gene_name = gene_node->getName();
      // species_candidates contains all candidates which are :
      // the index of the species and if the name is before the gene name

      for (auto &species_node: species_nodes) { // For each species_node (as std::shared_ptr<bpp::PhyloNode>)
        if(species_node->hasName()) { // if the node has a name
          std::string species_name = species_node->getName(); // get the name of the species
          if(species_name.size() <= gene_name.size()) { // if the gene name could contain the species name.
            if (Utils::string_comp(species_name, gene_name.substr(0, species_name.size()), case_sensitive_)) {
              // if the gene_name contains the species_name in the BEGINNING of the string,
              // push back in the species_candidates list */
              species_candidates.emplace_back(
                  std::make_pair(species_node,
                                 true) //true because of the position, before the gene name
              );
            }
            else if (Utils::string_comp(
                species_name
                , gene_name.substr(gene_name.size() - species_name.size(), gene_name.size())
                , case_sensitive_)) {
              // if the gene_name contains the species_name at the END of the string
              species_candidates.emplace_back( std::make_pair(species_node,
                                false) //false because of the position, after the gene name
                );
            }
          }
        }
      }

      // if there is multiple candidates, sort them by the name length
      // (longest first, smallest last)
      if(species_candidates.size() > 1) {
        species_candidates.sort([](const std::pair<std::shared_ptr<bpp::PhyloNode>, bool> &candidate1,
                                   const std::pair<std::shared_ptr<bpp::PhyloNode>, bool> &candidate2) {
          return candidate1.first->getName().size() > candidate2.first->getName().size();
        });
      }

      // There is some controls and warnings
      if (species_candidates.size() > 1) {
        std::cerr << "There are more than one species that could be associated "
                  << "with the gene, named: " << gene_name << std::endl;
        std::cerr << "\tCandidates are:" << std::endl;
        for (auto &pair: species_candidates) {
          std::cerr << "\t\t*" << pair.first->getName();
          if(pair.second){
            std::cerr << " placed at the beginning." << std::endl;
          } else {
            std::cerr << " placed at the end." << std::endl;
          }
        }
        std::cerr << "\tThe program has chosen "
                  << species_candidates.front().first->getName()
                  << " the species associated with the gene "
                  << gene_name << "." << std::endl << std::endl;
      }

      // Then, add the result to the map: a gene is associated to a species.
      if (not species_candidates.empty()) {
        genemap_returned.addGene(species_candidates.front().first, gene_node);
      } else {
        unmapped_gene_nodes.push_back(gene_node);
      }
    }
  }

  // Check if the unmapped gene nodes list is unmapped
  if(not unmapped_gene_nodes.empty()) {
    std::cerr << "Error during gene <> species mapping, "
              << "some gene leaves cannot be mapped:"
              << std::endl;
    Utils::write(std::cerr,
                 unmapped_gene_nodes.begin(), unmapped_gene_nodes.end(),
                 "", ", ", ".");
    std::cerr << std::endl;
    exit(EXIT_FAILURE);
  }

  return genemap_returned;
}

GeneMap<std::string, std::string>
SpeciesGeneMapper::mapWithTrees(const std::string& genetree, const std::string& speciestree) {

  // Initialization of the object returned by the function
  GeneMap<std::string, std::string> genemap_returned;

  // Get all nodes of each tree.
  auto species_tree = IO::newickToPhyloTree(speciestree);
  auto gene_tree = IO::newickToPhyloTree(genetree);

  auto genemap = mapWithTrees(*gene_tree, *species_tree);

  genemap_returned = nodeMapsToStrings(genemap, *gene_tree, *species_tree, true);

  return genemap_returned;
}

SpeciesGeneMap
SpeciesGeneMapper::mapFromFile(const std::string &filename, const bpp::PhyloTree &speciestree,
                               const bpp::PhyloTree &genetree) {
  SpeciesGeneMap ret = load(filename, speciestree.getAllLeaves(), genetree.getAllLeaves());
  //completeMap_(ret, speciestree, genetree);
  return ret;
}

GeneMap<std::string, std::string>
SpeciesGeneMapper::mapFromFile(const std::string& genetree
                               , const std::string& speciestree
                               , const std::string &filename
){
  auto gene_tree = IO::newickToPhyloTree(genetree);
  auto species_tree = IO::newickToPhyloTree(speciestree);

  auto map = mapFromFile(filename, *species_tree, *gene_tree);

  return nodeMapsToStrings(map, *gene_tree, *species_tree, true);
}

/* SAVE, LOAD FROM A FILE */

SpeciesGeneMap
SpeciesGeneMapper::load(const std::string &filename
                        , const std::vector<Node>& speciestree_nodes
                        , std::vector<Node> genetree_nodes
                        , bool printProgression) {

  // Open file
  std::ifstream file(filename, std::ios::in);

  // Test
  if(not file){
    std::cerr << "Error : map file " << filename << " does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Transform genetree_nodes into a umap with the node name as key.
  std::unordered_multimap<std::string, Node> gene_node_dictionary;
  for(auto gnode: genetree_nodes){
    std::string gname = gnode->getName();
    if(not case_sensitive_)
      std::transform(gname.begin(), gname.end(), gname.begin(), ::tolower);
    gene_node_dictionary.insert(std::pair<std::string, Node>(gname, gnode));
  }

  // Transform speciestree_nodes into a umap with the node name as key.
  std::unordered_map<std::string, Node> species_node_dictionary;
  for(auto snode: speciestree_nodes){
    std::string sname = snode->getName();
    if(not case_sensitive_)
      std::transform(sname.begin(), sname.end(), sname.begin(), ::tolower);
    species_node_dictionary[sname] = snode;
  }

  // Used to get line content
  std::string line_content;

  // Compute the number of lines.
  auto nlines = IO::nlines(file);

  std::size_t read_line = 0;
  SpeciesGeneMap ret;

  // For each line in file.
  while(std::getline(file, line_content) and not gene_node_dictionary.empty()) {

    #if defined _WIN32 || defined __CYGWIN__
    line_content.erase(
      std::remove(line_content.begin(), line_content.end(), '\r')
      , line_content.end());
    #endif

    // If not case_sensitive, change all characters for lowercase.
    if(not case_sensitive_)
      std::transform(line_content.begin(), line_content.end(), line_content.begin(), ::tolower);

    auto content = Utils::splitString(line_content, " \t");
    // If line content is not empty.
    if(not content.empty()) {
      // Species name is the second column.
      if(not (content.size() > 1)){
        std::cerr << "SMap file incorrect line " << read_line +1 << " read : " << line_content << std::endl;
      }
      std::string& species_name = content.at(1);


      // then associate species name to a species node from speciestree.
      Node species_node = nullptr;
      auto species_dict_it = species_node_dictionary.find(species_name);
      if(species_dict_it != species_node_dictionary.end()) {
        species_node = species_dict_it->second;
      }

      // If there is a species node associated with the species name, continue by looking for genes.
      if (species_node != nullptr) {

        // Get gene name pattern in the first column.
        std::string& gene_name_pattern = content.at(0);

        // Init gene_node which will be the associated gene node with the current species.
        Node gene_node = nullptr;

        // First, check if there is a perfect match with a gene in the gene_node_dictionary.
        auto gene_dict_it = gene_node_dictionary.find(gene_name_pattern);

        // If there is a perfect match, the gene name pattern is not a pattern but an exact name.
        if(gene_dict_it != gene_node_dictionary.end()) {
          // gene_name_pattern is an exact name.
          gene_node = gene_dict_it->second;
          ret.addGene(species_node, gene_node);
          gene_node_dictionary.erase(gene_dict_it);
        }
        // Else, it is a pattern, we will associate with the gene all names associated with the pattern.
        else {

          // We store all selected genes (matching with the pattern) into a list.
          std::list<decltype(gene_dict_it)> selected_genes;

          // For each gene node in the gene dictionnary, we will check if it matches with the pattern.

          for(gene_dict_it = gene_node_dictionary.begin(); gene_dict_it != gene_node_dictionary.end(); gene_dict_it++){
            auto& gene_name = gene_dict_it->first; // get gene name (the key of our dictionary).
            auto& gene_node = gene_dict_it->second; // get gene node in tree (the value of our dictionnary).

            // If the gene matches with the pattern, associate with the species.
            if(PhyloTreeToolBox::node_name_match_with_pattern(gene_name, gene_name_pattern, case_sensitive_)){
              ret.addGene(species_node, gene_node);
              selected_genes.push_back(gene_dict_it);
            } else {}
          }
          // then delete all associated genes from dictionary.
          for(auto& it: selected_genes) gene_node_dictionary.erase(it);
        }
      } else
        std::cerr << __FUNCTION__ << " warning: " << species_name << " does not exist in the speciestree." << std::endl;
    }
    read_line++;
  }

  // Everything is ok if all gene is associated with a species -> the dictionary is empty.
  if(not gene_node_dictionary.empty()) {
    std::cerr << "Error during gene <> species mapping, some gene leaves cannot be mapped:" << std::endl;
    for(const auto& unmapped_gene: gene_node_dictionary) {
      std::cout << "* " << unmapped_gene.second << std::endl;
    }
    std::cerr << "Please check if the map file " << filename << " is correct." << std::endl;
    exit(EXIT_FAILURE);
  }
  assert(gene_node_dictionary.empty());

  file.close();
  return ret;

}

SpeciesGeneMap
SpeciesGeneMapper::loadCompara(const std::string &filename
                               , const std::vector<Node>& speciestree_nodes
                               , std::vector<Node> genetree_nodes
                               , const bool printProgression) {
  // Open file
  std::ifstream file(filename, std::ios::in);

  // Test
  if(not file){
    std::cerr << "Error : map file " << filename << " does not exist." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Transform genetree_nodes into a mmap with the node name as key.
  std::unordered_multimap<std::string, Node> gene_node_dictionnary;
  for(auto gnode: genetree_nodes){
    std::string gname = gnode->getName();
    if(not case_sensitive_)
      std::transform(gname.begin(), gname.end(), gname.begin(), ::tolower);
    gene_node_dictionnary.insert(std::pair<std::string, Node>(gname, gnode));
  }

  // Transform speciestree_nodes into a umap with the node name as key.
  std::unordered_map<std::string, Node> species_node_dictionnary;
  for(auto snode: speciestree_nodes){
    std::string sname = snode->getName();
    if(not case_sensitive_)
      std::transform(sname.begin(), sname.end(), sname.begin(), ::tolower);
    species_node_dictionnary[sname] = snode;
  }

  // Used to get line content
  std::string line_content;

  // Compute the number of lines
  auto nlines = IO::nlines(file);

  std::size_t read_line = 0;
  SpeciesGeneMap ret;
  while(std::getline(file, line_content) and not gene_node_dictionnary.empty()){
    #if defined _WIN32 || defined __CYGWIN__
    file_content.erase(
      std::remove(line_content.begin(), line_content.end(), '\r')
      , line_content.end());
    #endif

    auto content = Utils::splitString(line_content, " ");
    if(not content.empty()) {
      if (content.at(0).compare("SEQ") == 0) { // If the line is a mapping line.
        std::string& species_name = content.at(1);
        if(not case_sensitive_)
          std::transform(species_name.begin(), species_name.end(), species_name.begin(), ::tolower);


        // then associate species name to a species node from speciestree.
        Node species_node = nullptr;
        auto species_dict_it = species_node_dictionnary.find(species_name);
        if(species_dict_it != species_node_dictionnary.end()) {
          species_node = species_dict_it->second;
        }

        if (species_node != nullptr) {
          std::string& gene_name = content.at(2);
          if(not case_sensitive_)
            std::transform(gene_name.begin(), gene_name.end(), gene_name.begin(), ::tolower);
          Node  gene_node = nullptr;
          auto gene_dict_it = gene_node_dictionnary.find(gene_name);
          if(gene_dict_it != gene_node_dictionnary.end()) {
            gene_node = gene_dict_it->second;
            ret.addGene(species_node, gene_node);
            gene_node_dictionnary.erase(gene_dict_it);
          } else { }
        } else
          std::cerr << __FUNCTION__ << " error: " << species_name << " does not exist in the speciestree." << std::endl;
      }
    }
    read_line++;
  }

  std::list<Node> unmapped_nodes {};
  for (auto& unmapped_node: gene_node_dictionnary) {
    unmapped_nodes.push_back(unmapped_node.second);
  }

  if (not unmapped_nodes.empty()) {
    std::cerr << "Error during gene <> species mapping, "
              << "some gene leaves cannot be mapped:" << std::endl;
    Utils::write(std::cerr, unmapped_nodes.begin(), unmapped_nodes.end(),
                 "", ", ", ".");
    std::cerr << std::endl;
    exit(EXIT_FAILURE);
  }

  file.close();
  return ret;
}

GeneMap<std::size_t, std::size_t>
SpeciesGeneMapper::nodeMapsToIndexes(const SpeciesGeneMap &genemap,
                                     const bpp::PhyloTree &genetree,
                                     const bpp::PhyloTree &speciestree) {
  GeneMap<std::size_t, std::size_t> ret;

  auto map = genemap.u_map();

  for(auto it = map.begin(); it != map.end(); it++){
    std::vector<std::size_t> genes_indexes;
    for(auto& gene: it->second){
      genes_indexes.emplace_back(genetree.getNodeIndex(gene));
    }
    ret.addGenes(speciestree.getNodeIndex(it->first), genes_indexes);
  }

  return ret;
}

GeneMap<std::string, std::string> SpeciesGeneMapper::nodeMapsToStrings(const SpeciesGeneMap &genemap) {
  GeneMap<std::string, std::string> ret;

  // Get all species nodes.
  auto species_nodes = genemap.getSpecies();

  // Then get all genes associated with these species.
  // A node (gene as species) will be inserted in the resulting genemap if it has a name.
  for(auto& species: species_nodes){
    if(species->hasName()){
      auto associated_genes = genemap.getGenes(species);
      for(auto& gene: associated_genes){
        if(gene->hasName()){
          ret.addGene(species->getName(), gene->getName());
        }
      }
    }
  }

  return ret;
}

GeneMap<std::string, std::string> SpeciesGeneMapper::nodeMapsToStrings(const SpeciesGeneMap &genemap
                                                  , const bpp::PhyloTree &genetree
                                                  , const bpp::PhyloTree &speciestree
                                                  , const bool leaves_only) {
  GeneMap<std::string, std::string> ret;

  // Get all species nodes.
  std::list<decltype(speciestree.getRoot())> species_nodes;
  if(leaves_only)
    species_nodes = PhyloTreeToolBox::getLeaves(speciestree);
  else {
    species_nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalRecursive_list(speciestree);
  }

  // Then get all genes associated with these species.
  // A node (from genetree as speciestree) will be inserted in the resulting genemap if it has a name and if it follows
  // the condition leaves_only.
  for(auto& species: species_nodes){
    if(species->hasName()){
      auto associated_genes = genemap.getGenes(species);
      for(auto& gene: associated_genes){
        if(gene->hasName() and not (leaves_only and not genetree.isLeaf(gene))){
          // not (leaves_only and not genetree.isLeaf(gene)) is NAND gate
          // ____________________________________________________
          // | leaves_only | not genetree.isLeaf(gene) | result |
          // ----------------------------------------------------
          // |       false |                    false  |   true |
          // |       false |                     true  |   true |
          // |        true |                    false  |   true |
          // |        true |                     true  |  false |
          // ----------------------------------------------------
          ret.addGene(species->getName(), gene->getName());
        }
      }
    }
  }

  return ret;
}

SpeciesGeneMap SpeciesGeneMapper::reduceMapAccordingToGeneList(
    const SpeciesGeneMap &map,
    const std::vector<std::shared_ptr<bpp::PhyloNode>>& genes) {
  SpeciesGeneMap result;
  for(auto& species: map.getSpecies()){
    for(auto it_gene = genes.begin(); it_gene != genes.end(); it_gene++){
      if(map.speciesHasGene(species, *it_gene)){
        result.addGene(species, *it_gene);
      }
    }
  }
  return result;
}

SpeciesGeneMap SpeciesGeneMapper::aggregate(
    const SpeciesGeneMap &map1, const SpeciesGeneMap &map2) {
  SpeciesGeneMap res(map1);
  res.addMap(map2);
  return res;
}

void SpeciesGeneMapper::updateInnerNodeMapping(
    SpeciesGeneMap &genemap, const bpp::PhyloTree &genetree,
    const bpp::PhyloTree &speciestree) {
  auto genetree_inner_nodes = PhyloTreeToolBox::getInternalNodes(genetree);

  //auto genetree_inner_nodes = genetree.getAllInnerNodes();
  for(auto& ginode: genetree_inner_nodes){
    associate_species_with_ancestral_gene_(genemap, ginode,
                                           genetree, speciestree);
  }
}

void SpeciesGeneMapper::updateInnerNodeMappingAfterRerooting(
    SpeciesGeneMap& genemap,
    const bpp::PhyloTree& genetree, const bpp::PhyloTree& speciestree,
    const Node& old_root) {
  auto current_node = old_root;
  std::list<decltype(current_node)> path_to_root {old_root};
  while(genetree.hasFather(current_node)) {
    current_node = genetree.getFather(current_node);
    path_to_root.push_back(current_node);
  }

  for(auto& ginode: path_to_root){
    associate_species_with_ancestral_gene_(genemap, ginode,
                                           genetree, speciestree);
  }
};

bool SpeciesGeneMapper::checkIfAllGenesAreMapped(
    const SpeciesGeneMap &map, const bpp::PhyloTree &genetree
    , const bpp::PhyloTree &speciestree, bool gene_leaves_only, bool quiet
) {
  auto gene_nodes = (gene_leaves_only) ? genetree.getAllLeaves() : genetree.getAllNodes();
  for(auto gene_node: gene_nodes) {
    if(not map.hasGene(gene_node)){
      if(not quiet) {
        std::cerr << "Error in species <> gene mapping, incorrect data."
                  << std::endl;
        std::cerr << "Genes are missing in map." << std::endl;
        std::cerr << "Current map:" << std::endl << map.u_map() << std::endl;
        std::cerr << "Missing gene tree nodes:" << std::endl;
        for (auto gene_node_: gene_nodes) {
          if (not map.hasGene(gene_node_)) {
            std::cerr << "\t* " << gene_node_ << std::endl;
          }
        }
      }
      return false;
    } else {
      if(not map.getAssociatedSpecies(gene_node)) {
        if(not quiet) {
          std::cerr << "Error in species <> gene mapping, incorrect data."
                    << std::endl;
          std::cerr << "Some genes are not associated in map with a species."
                    << std::endl;
          std::cerr << "Current map:" << std::endl << map.u_map() << std::endl;
          std::cerr << "Related to genes:" << std::endl;
          for (auto gene_node_: gene_nodes) {
            if (not map.getAssociatedSpecies(gene_node_)) {
              std::cerr << "\t* " << gene_node_ << std::endl;
            }
          }
        }
        return false;
      }
    }

    if(not map.speciesHasGene(map.getAssociatedSpecies(gene_node), gene_node)){
      if(not quiet) {
        std::cerr << "Error in species <> gene mapping, incorrect data."
                  << std::endl;
        std::cerr << "Error during map building." << std::endl;
        std::cerr << "Current associated species with " << gene_node << ": "
                  << map.getAssociatedSpecies(gene_node)
                  << std::endl;
        std::cerr << "Current map:" << std::endl << map.u_map() << std::endl;
      }
      return false;
    }
  }
  return true;
}

void
SpeciesGeneMapper::save(const std::string &filename, const SpeciesGeneMap &genemap, const bpp::PhyloTree &speciestree,
                        const bpp::PhyloTree &genetree) {
  std::ofstream file(filename);

  writeSMap(file, genemap, speciestree, genetree, false, false);

  file.close();
}

void SpeciesGeneMapper::save(const std::string &filename
                             , const std::vector<SpeciesGeneMap> &genemaps
                             , const bpp::PhyloTree &speciestree
                             , const std::vector<std::shared_ptr<bpp::PhyloTree>> &genetrees
                             , const bool printProgression)
{
  assert(genemaps.size() == genetrees.size());

  std::ofstream file(filename);


  for(std::size_t i = 0; i < genetrees.size(); ++i){
    auto genetree = genetrees.at(i);

    const SpeciesGeneMap& genemap = genemaps.at(i);

    writeSMap(file, genemap, speciestree, *genetree, false, true);

  }

  file.close();
}

void SpeciesGeneMapper::writeSMap(
    std::ostream &os
    , const SpeciesGeneMap &map
    , const bpp::PhyloTree &speciestree
    , const bpp::PhyloTree &genetree
    , const bool printProgression
    , const bool writeInternalSpecies
) {
  auto genetree_leaves = genetree.getAllLeaves();


  for(const auto& genetree_leaf: genetree_leaves){
    if(map.hasGene(genetree_leaf)){
      auto associated_species = map.getAssociatedSpecies(genetree_leaf);
      if(speciestree.isLeaf(associated_species) or writeInternalSpecies) {
        assert(genetree_leaf->hasName());
        assert(associated_species->hasName());
        os << genetree_leaf->getName() << " " << associated_species->getName() << std::endl;
      }
    } else {
      std::cerr << "Warning: gene " << genetree_leaf << " is not mapped" << std::endl;
      assert(not map.hasGene(genetree_leaf));
    }
  }

}

void SpeciesGeneMapper::associate_species_with_ancestral_gene_(
    SpeciesGeneMap &map, const Node &node, const bpp::PhyloTree &genetree
    , const bpp::PhyloTree &speciestree
) {
  auto gene_sons = genetree.getSons(node);

  // Find species associated with the sons of the inner node.
  std::vector<Node> species_sons(gene_sons.size());
  std::size_t i = 0;
  std::generate(species_sons.begin(), species_sons.end(),
                [&i, &map, &gene_sons]{
                  return map.getAssociatedSpecies(gene_sons.at(i++)); }
  );

  // The species associated to the current inner_node is the common species ancestor of these species.
  Node last_common_ancestor
      = PhyloTreeToolBox::getLastCommonAncestor(species_sons, speciestree);
  map.setGene(last_common_ancestor, node);
}
