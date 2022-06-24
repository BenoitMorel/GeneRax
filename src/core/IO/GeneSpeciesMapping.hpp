#pragma once

#include <string>
#include <map>
#include <fstream>
#include <unordered_set>

typedef struct corax_utree_s corax_utree_t;
typedef struct corax_unode_s corax_unode_t;
typedef struct corax_rtree_s corax_rtree_t;



/**
 *  Parse and give access to gene-to-species mapping
 *  Support the following format:
 *    species1:gene1;gene2;gene3
 *    species2:gene1
 *    species3:gene1:gene2;gene3
 *
 */
class GeneSpeciesMapping {
public:
  
  GeneSpeciesMapping() {}
  
  void fill(const std::string &mappingFile, const std::string &geneTreeStrOrFile); 
  void fillFromGeneLabels(const std::unordered_set<std::string> &labels); 
  bool check(corax_utree_t *geneTree, corax_rtree_t *speciesTree);
  bool check(const std::unordered_set<std::string> &genes, const std::unordered_set<std::string> &species);
  void fill(const GeneSpeciesMapping &mapping);

  /**
   *  @return an object mapping each gene to one species
   */
  const std::map<std::string, std::string> &getMap() const {return _map;}

  /**
   *  @param gene gene std::string
   *  @return the species mapped to this gene
   */
  const std::string &getSpecies(const std::string &gene) const {return _map.find(gene)->second;}

  std::unordered_set<std::string> getCoveredSpecies() const;
private:
  std::map<std::string, std::string> _map; // <gene,species>
  void buildFromPhyldogMapping(std::ifstream &f);
  void buildFromTreerecsMapping(std::ifstream &f);
  void buildFromMappingFile(const std::string &mappingFile); 
  void buildFromTrees(const std::string &geneTreeStrOrFile); 
};

