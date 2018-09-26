// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 14/12/16.
//

#ifndef PHYLASOLVER_SPECIESGENEMAPPER_H
#define PHYLASOLVER_SPECIESGENEMAPPER_H

// cpp libs
#include <Bpp/Phyl/Io/Newick.h>
#include <cassert>
#include <string>

// Treerecs libs
#include "PhyloTreeToolBox.h"
#include <ale/containers/GeneMap.h>


/// Defines the mapping method.
enum MAPPING_METHOD {
  trees /// the automatic way using only trees (default).
  , nhxtag /// use NHX tags.
  , genenames /// according to the gene names.
  , smap /// with a smap file.
  , ensembl /// with an ensembl file.
};

/*!
 * \class SpeciesGeneMapper
 * \brief SpeciesGeneMapper provides methods to build a map according to a mapping file, or by the names in gene and
 *        species tree. It also provides methods to edit or translate a SpeciesGeneMap.
 * \details SpeciesGeneMapper contains 3 mapping methods :
 *          * According to gene names by specification of a separator (with position) between species and gene name in gene tree leaf names: `mapWithGeneNames(genetree, speciestree, separator, species_name_is_before_separator)` ;
 *          * According to a file by specifying a smap or emf file: `mapFromFile(filename, speciestree, genetree)` ;
 *          * According to the trees only: SpeciesGeneMapper will seek in gene names species names according to the species tree leaves `mapWithTrees(genetree, speciestree)`.
 *
 *          SpeciesGeneMapper can also predict internal genes association with species nodes according to the last species common ancestor.
 */

class SpeciesGeneMapper {
protected:

  /// Define for each comparison of characters if it's case sensitive or not. See Constants::DEFAULT_CASE_SENSITIVE_MAPPING.
  static bool case_sensitive_;

  /// \brief Extracts the species name from the gene name according to a separator the position before (true) or after (false)
  ///        the separator.
  /// \param before Position of the species name before the separator.
  /// \return The species name associated with the gene.
  static std::string get_specie_from_gene_name_(const std::string &genename, const std::string& separator = "_", bool before = true) ;

  /// \brief Update an ancestral node map according to the species affiliation of its sons.
  static void associate_species_with_ancestral_gene_(SpeciesGeneMap& map, const Node& node, const bpp::PhyloTree& genetree, const bpp::PhyloTree& speciestree);

public:
/****************
 * Constructor
 */
  //SpeciesGeneMapper(const bool case_sensitive = false): case_sensitive_(case_sensitive) {};

/****************
 * Destructor
 */
  ~SpeciesGeneMapper() = default;

/****************
 * Methods
 */


  /*!
   * @brief Create Species <> gene maps according to a MAPPING_METHOD.
   * @details There is three methods to generate a SpeciesGeneMap:
   *        * Reading a file which contains correspondences (SMap format usually);
   *        * Looking in gene names for the species names (with a separator and a given position);
   *        * Idem, but automatically looking the names in species tree and gene tree (default).
   *
   *        See MAPPING_METHOD enum.
   */
  template <typename GeneTreeIterator>
  static std::vector<SpeciesGeneMap> map(
      const GeneTreeIterator& genetrees_begin /// Iterator to the first genetree to use.
      , const GeneTreeIterator& genetrees_end /// Iterator to the last genetree to use.
      , const bpp::PhyloTree& speciestree /// Species tree.
      , const MAPPING_METHOD& method /// Mapping method.
      , const std::string& mapfilename = "" /// Map filename if the method uses a file.
      , const std::string& separator = DEFAULT_SEPARATION_CHARACTER_BETWEEN_GENE_AND_SPECIES_NAME /// Character which separates species name to the gene name.
      , const bool& species_position_is_prefix = DEFAULT_SPECIES_NAME_IS_PREFIX_IN_GENE_NAME /// Position (prefix or not) of the species name according to the separator.
      , const bool printProgression = true /// Print progression bar.
  ) {

    std::vector<std::shared_ptr<bpp::PhyloTree>> genetrees {genetrees_begin, genetrees_end};
    auto genemaps = std::vector<SpeciesGeneMap>(genetrees.size());
    if(method == ensembl or method == smap){
      // First collect in nodes all leaves of all genetrees.
      std::vector<Node> nodes;
      for(auto& genetree: genetrees){
        std::vector<Node> genetree_nodes = genetree->getAllLeaves();
        nodes.insert(nodes.end(), genetree_nodes.begin(), genetree_nodes.end());
      }

      // Then map all of theses in a global map.
      SpeciesGeneMap global_genemap;
      if (method == ensembl) {
        global_genemap = SpeciesGeneMapper::loadCompara(mapfilename, speciestree.getAllLeaves(), nodes, printProgression);
      } else if (method == smap) {
        global_genemap = SpeciesGeneMapper::load(mapfilename, speciestree.getAllLeaves(), nodes, printProgression);
      }

      // Split global map into maps for each gene tree increases performances.
      genemaps = SpeciesGeneMapper::splitGlobalMap(genetrees.begin(), genetrees.end(), global_genemap);

      // End of mapping with files.

    } else {
      // Mapping with gene names.
      // * automatically by finding the species name in the gene name.
      // * manually with specification of species name position (prefix/postfix) and separator.
      #if defined(_OPENMP)
      # pragma omp parallel for schedule(dynamic)
      #endif
      for (std::size_t i = 0; i < genetrees.size(); i++) {
        auto& genetree = genetrees.at(i);
        SpeciesGeneMap genetree_genemap;
        if (method == genenames) {
          // Using a separator and species name position (prefix/postfix)
          genetree_genemap
              = SpeciesGeneMapper::mapWithGenesNames(
              *genetree, speciestree,
              separator, species_position_is_prefix);
        } else if (method == nhxtag) {
          // Use NHX tags informations (if gene tree has been written in NHX
          // with "S" tag).
          genetree_genemap
              = SpeciesGeneMapper::mapWithNhxTags(*genetree, speciestree);
        } else {
          // Auto
          assert(method == trees);
          genetree_genemap
              = SpeciesGeneMapper::mapWithTrees(*genetree, speciestree);
        }
        #if defined(_OPENMP)
        # pragma omp critical
        #endif
        {
          genemaps[i].addMap(genetree_genemap);
        }
      }
      // End of mapping with gene names.
    }

    // Test if every map is correct. Raise a warning or an error when there is something wrong.
    if(not SpeciesGeneMapper::checkSpeciesGeneMapsConsistence(genemaps.begin(), genemaps.end(), genetrees.begin(), genetrees.end(), speciestree)) {
      exit(EXIT_FAILURE);
    }

    return genemaps;
  }

  /// \brief Check if genes and species can be associated according to nhx tags.
  /// \return True if gene tree can be associated with the species tree
  ///         according to nhx tags.
  static bool mappableWithNhxTags(const bpp::PhyloTree &genetree,
      const bpp::PhyloTree &speciestree,
      const std::string &tag = "Species name");

  /// \brief Check if genes and species can be associated according to nhx tags.
  /// \return True if gene tree can be associated with the species tree
  ///         according to nhx tags.
  template<typename GeneTreesIterator>
  static bool mappableWithNhxTags(
      GeneTreesIterator genetrees_begin,
      const GeneTreesIterator& genetrees_end,
      const bpp::PhyloTree &speciestree,
      const std::string &tag = "Species name") {
    for(;genetrees_begin != genetrees_end; genetrees_begin++) {
      if(not mappableWithNhxTags(**genetrees_begin, speciestree, tag))
        return false;
    }
    return true;
  };

  /* MAPPING with a map of nodes returned */

  /// \brief Species to genes mapping according to NHX tags in gene trees
  /// \param tag Name of the tag to get species name in gene node.
  /// \return A Genemap<std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>> which contains relations
  ///         between species and genes.
  static SpeciesGeneMap mapWithNhxTags(
      const bpp::PhyloTree &genetree,
      const bpp::PhyloTree &speciestree,
      const std::string &tag = "Species name");

  /// \brief Species to genes mapping according to NHX tags in gene trees
  /// \param tag Name of the tag to get species name in gene node.
  /// \return A Genemap<std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>> which contains relations
  ///         between species and genes..
  static GeneMap<std::string, std::string> mapWithNhxTags(
      const std::string &genetree,
      const std::string &speciestree,
      const std::string &tag = "Species name");

  /// \brief Species to genes mapping according to the species names extracted from the gene names.
  /// \param before Seek species name before the separator.
  /// \return A Genemap<std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>> which contains relations
  ///         between species and genes.
  static SpeciesGeneMap mapWithGenesNames(
      const bpp::PhyloTree &genetree,
      const bpp::PhyloTree &speciestree,
      std::string separator = "_",
      const bool before = true);

  /// \brief Species to genes mapping according to the species names extracted from the gene names.
  ///        If before = true (default), species name is before the separator.
  /// \param before Seek species name before the separator.
  /// \return A Genemap<std::string, std::string> which contains relations
  ///         between species and gene leaf names.
  static GeneMap<std::string, std::string> mapWithGenesNames(
      const std::string &genetree,
      const std::string &speciestree,
      std::string separator = "_",
      const bool before = true /// Seek species name before the separator.
  );

  /// \brief Species to genes mapping according to the species tree's leaves names and the gene tree's leaves names. The names
  ///        extraction is automatic (no separator or position given).
  /// \return A Genemap<std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>> which contains relations
  ///         between species and genes.
  static SpeciesGeneMap mapWithTrees(
      const bpp::PhyloTree &genetree, const bpp::PhyloTree &speciestree);

  /// \brief Species to genes mapping according to the species tree's leaves names and the gene tree's leaves names. The names
  /// extraction is automatic (no separator or position given).
  /// \return A Genemap<std::string, std::string> which contains relations
  ///         between species and gene leaf names.
  static GeneMap<std::string, std::string>
  mapWithTrees(const std::string& genetree, const std::string& speciestree);

  /// \brief Read a SpeciesGeneMap file and complete it according to the speciestree and genetree.
  /// \return A Genemap<std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>> which contains relations
  ///         between species and genes.
  static SpeciesGeneMap mapFromFile(
      const std::string &filename,
      const bpp::PhyloTree &speciestree,
      const bpp::PhyloTree &genetree);

  /// \brief Read a SpeciesGeneMap file and complete it according to the speciestree and genetree std::strings.
  /// \return A Genemap<std::string, std::string> which contains relations
  ///         between species and gene leaf names.
  static GeneMap<std::string, std::string> mapFromFile(
      const std::string& genetree
    , const std::string& speciestree
    , const std::string& filename);

  /// \brief Creates a SpeciesGeneMap of indexes (std::size_t) instead of a map of std::shared_ptr<bpp::PhyloNode>.
  /// \return A Genemap<std::size_t, std::size_t> which contains relations
  ///         between species and gene leaf indexes in each tree.
  static GeneMap<std::size_t, std::size_t> nodeMapsToIndexes(const SpeciesGeneMap &genemap,
                                                      const bpp::PhyloTree &genetree,
                                                      const bpp::PhyloTree &speciestree);

  /// \brief Creates a SpeciesGeneMap of names (std::string) instead of a map of std::shared_ptr<bpp::PhyloNode>.
  /// \return A Genemap of names (including internal node names if the given SpeciesGeneMap contains that).
  static GeneMap<std::string, std::string> nodeMapsToStrings(const SpeciesGeneMap &genemap);

  /// \brief Creates a SpeciesGeneMap of names (std::string) instead of a map of std::shared_ptr<bpp::PhyloNode>.
  /// \param leaves_only keep only leaf names in the resulting map.
  /// \return A Genemap of names.
  static GeneMap<std::string, std::string> nodeMapsToStrings(const SpeciesGeneMap &genemap,
                                                      const bpp::PhyloTree &genetree,
                                                      const bpp::PhyloTree &speciestree,
                                                      const bool leaves_only = true);

  /* LOAD, SAVE a map to a file */

  /// \brief Read smap file and returns a map of nodes.
  static SpeciesGeneMap load(const std::string& filename, const std::vector<Node>& speciestree_nodes, std::vector<Node> genetree_nodes, bool printProgression = true);

  /// \brief Read a Compara file and returns a map of nodes.
  static SpeciesGeneMap loadCompara(const std::string& filename, const std::vector<Node>& speciestree_nodes, std::vector<Node> genetree_nodes, bool printProgression = true) ;

  /// \brief Save a map of a given genetree.
  static void writeSMap(std::ostream& os
      , const SpeciesGeneMap& map
      , const bpp::PhyloTree& speciestree
      , const bpp::PhyloTree& genetree
      , const bool printProgression = true
      , const bool writeInternalSpecies = true
  );

  /// \brief Write a map file in the smap standard.
  static void save(const std::string& filename
                   , const SpeciesGeneMap& genemap
                   , const bpp::PhyloTree& speciestree
                   , const bpp::PhyloTree& genetree);

  /// \brief Write a map file in the smap standard.
  template<typename Genetrees_iterator>
  static void save(const std::string &filename,
                          const SpeciesGeneMap &genemap,
                          const bpp::PhyloTree& speciestree,
                          const Genetrees_iterator& genetrees_begin,
                          const Genetrees_iterator& genetrees_end,
                          const bool printProgression = true
  ) {
    std::ofstream file(filename);


    for(auto genetree_iterator = genetrees_begin; genetree_iterator != genetrees_end; genetree_iterator++){
      const auto& genetree = *genetree_iterator;
      writeSMap(file, genemap, speciestree, *genetree, false, true);
    }

    file.close();
  }

  /// \brief Write a vector of maps in a file in the smap standard.
  static void save(const std::string &filename,
            const std::vector<SpeciesGeneMap>& genemaps,
            const bpp::PhyloTree& speciestree,
            const std::vector<std::shared_ptr<bpp::PhyloTree>>& genetrees,
            const bool printProgression = true
  );

  /// \brief Split a map (which contains mapping of a multiple gene trees) into a collection of maps (std::vector).
  /// \return A std::vector of SpeciesGeneMaps for each corresponding gene tree.
  template<typename Genetree_iterator>
  static std::vector<SpeciesGeneMap> splitGlobalMap(
      const Genetree_iterator& gtrees_begin
      , const Genetree_iterator gtrees_end
      , const SpeciesGeneMap& global_map
  ) {
    assert(std::distance(gtrees_begin, gtrees_end) != 0);

    std::vector<SpeciesGeneMap> genemaps;
    genemaps.reserve((unsigned long) std::distance(gtrees_begin, gtrees_end));

    for(auto gtree_it = gtrees_begin; gtree_it != gtrees_end; gtree_it++){
      SpeciesGeneMap genemap;

      std::shared_ptr<bpp::PhyloTree> gtree = *gtree_it;
      auto gtree_leaves = gtree->getAllLeaves();

      for(auto& gtree_leaf: gtree_leaves){
        auto gtree_leaf_species = global_map.getAssociatedSpecies(gtree_leaf);
        genemap.addGene(gtree_leaf_species, gtree_leaf);
      }

      genemaps.emplace_back(std::move(genemap));
    }

    return genemaps;
  }

  /// \brief Reduce Map according to a vector of gene nodes.
  /// \return A new SpeciesGeneMap which contains nodes of the genes vector only.
  static SpeciesGeneMap reduceMapAccordingToGeneList(
      const SpeciesGeneMap& map
      , const std::vector<std::shared_ptr<bpp::PhyloNode>>& genes
  );

  /// \brief Aggregate two maps.
  /// \return A new SpeciesGeneMap which contains informations of the two given maps.
  static SpeciesGeneMap aggregate(const SpeciesGeneMap& map1, const SpeciesGeneMap& map2) ;

  /// \brief Update a map by association of species nodes with internal gene names according to the las common ancestor.
  /// \return Nothing by modify the given genemap.
  static void updateInnerNodeMapping(SpeciesGeneMap& genemap, const bpp::PhyloTree& genetree, const bpp::PhyloTree& speciestree);

  /// \brief Update a map after a rerooting.
  /// \return Nothing by modify the given genemap.
  static void updateInnerNodeMappingAfterRerooting(SpeciesGeneMap& genemap, const bpp::PhyloTree& genetree, const bpp::PhyloTree& speciestree, const Node& old_root);

  /// \brief Set mapping case sensitive.
  /// \return Nothing but set next mappings case sensitive.
  static void set_case_sensitive_mapping() {
    SpeciesGeneMapper::case_sensitive_ = true;
  }

  /// \brief Set mapping non sensitive.
  /// \return Nothing but set next mappings case non sensitive.
  static void set_case_insensitive_mapping() {
    SpeciesGeneMapper::case_sensitive_ = false;
  }

  /// \brief Check if mapping is currently case sensitive.
  static bool case_sensitive() {
    return SpeciesGeneMapper::case_sensitive_;
  }

  /// \brief Check SpeciesGeneMap consistence with trees.
  static bool checkIfAllGenesAreMapped(
      const SpeciesGeneMap &map /// map to test.
      , const bpp::PhyloTree &genetree /// genetree associated with the map.
      , const bpp::PhyloTree &speciestree /// speciestree associated with the map.
      , bool gene_leaves_only = true /// check only with gene tree leaves.
      , bool quiet = true /// silent mode.
  );

  /// \brief Check SpeciesGeneMaps consistence with trees.
  /// \details Check if SpeciesgeneMaps are corrects with genetrees and speciestrees.
  ///          Note: returns false also in the case where the number of elements is not respected (number of maps == number of genetrees for example).
  template <typename MapsIterator, typename GenetreesIterator>
  static bool checkSpeciesGeneMapsConsistence(
      MapsIterator maps_begin /// begin of the map to test.
      , const MapsIterator& maps_end /// end of the map to test.
      , GenetreesIterator genetrees_begin /// begin of gene trees to test.
      , const GenetreesIterator& genetrees_end /// end of gene trees to test.
      , const bpp::PhyloTree& speciestree /// species tree.
      , bool gene_leaves_only = true /// check only with gene tree leaves.
  ){
    auto n_maps = std::distance(maps_begin, maps_end);
    auto n_genetrees = std::distance(genetrees_begin, genetrees_end);

    if(n_maps != n_genetrees){
      return false;
    }

    while(maps_begin != maps_end) {
      bool result_for_one_map = SpeciesGeneMapper::checkIfAllGenesAreMapped(
          *maps_begin, **genetrees_begin, speciestree, gene_leaves_only);

      if(not result_for_one_map) return false;

      // Get next elements
      maps_begin++;
      genetrees_begin++;
    }

    return true;
  }
};

#endif //PHYLASOLVER_SPECIESGENEMAPPER_H
