#include "GeneSpeciesMapping.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>


GeneSpeciesMapping::GeneSpeciesMapping(const string &mappingFile)
{
  ifstream f(mappingFile);
  string line;
  getline(f, line);
  f.seekg(0, ios::beg);
  if (line.find(":") != string::npos) {
    buildFromPhyldogMapping(f);
  } else {
    buildFromTreerecsMapping(f);
  }
}


void GeneSpeciesMapping::buildFromPhyldogMapping(ifstream &f)
{
  /*
    species1:gene1;gene2;gene3
    species2:gene4;gene5
  */
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

void GeneSpeciesMapping::buildFromTreerecsMapping(ifstream &f)
{
  /*
    gene1 species1
    gene2 species2
    gene3 species1
  */
  string line;
  while (getline(f, line)) {
    stringstream ss(line);
    string species;
    string gene;
    ss >> gene;
    ss >> species;
    _map[gene] = species;
  }
}

