#pragma once

#include <string>
#include <map>
#include <fstream>



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
  
  void fill(const std::string &mappingFile, const std::string &geneTreeStr); 

  /**
   *  @return an object mapping each gene to one species
   */
  const std::map<std::string, std::string> &getMap() const {return _map;}

  /**
   *  @param gene gene std::string
   *  @return the species mapped to this gene
   */
  const std::string &getSpecies(const std::string &gene) {return _map[gene];}

private:
  std::map<std::string, std::string> _map;
  void buildFromPhyldogMapping(std::ifstream &f);
  void buildFromTreerecsMapping(std::ifstream &f);
  void buildFromMappingFile(const std::string &mappingFile); 
  void buildFromTrees(const std::string &geneTreeStr); 
};

