#pragma once

#include <string>
#include <map>
#include <fstream>

using namespace std;

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
  
  void fill(const string &mappingFile, const string &geneTreeStr); 

  /**
   *  @return an object mapping each gene to one species
   */
  const map<string, string> &getMap() const {return _map;}

  /**
   *  @param gene: gene string
   *  @return the species mapped to this gene
   */
  const string &getSpecies(const string &gene) {return _map[gene];}

private:
  map<string, string> _map;
  void buildFromPhyldogMapping(ifstream &f);
  void buildFromTreerecsMapping(ifstream &f);
  void buildFromMappingFile(const string &mappingFile); 
  void buildFromTrees(const string &geneTreeStr); 
};

