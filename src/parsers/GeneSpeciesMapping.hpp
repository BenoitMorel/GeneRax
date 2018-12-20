#pragma once

#include <string>
#include <map>

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
  /**
   *  Constructor
   *  @param mappingFile: file containing the mapping information
   */
  GeneSpeciesMapping(const string &mappingFile);

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
};

