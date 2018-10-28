#include "GeneSpeciesMapping.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

GeneSpeciesMapping::GeneSpeciesMapping(const string &mappingFile)
{
  ifstream f(mappingFile);
  string line;
  while (getline(f, line)) {
    stringstream ss(line);
    string species;
    string gene;
    getline(ss, species, ':');
    while(getline(ss, gene, ';')) {
      _map[gene] = species;
    }
  }
}
