#include "GeneSpeciesMapping.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>

void GeneSpeciesMapping::fill(const string &mappingFile, const string &geneTreeStr) 
{
  if (mappingFile.size()) {
    buildFromMappingFile(mappingFile);
  } else {
    buildFromTrees(geneTreeStr);
  }
}

void GeneSpeciesMapping::buildFromMappingFile(const string &mappingFile)
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

void GeneSpeciesMapping::buildFromTrees(const string &geneTreeStr)
{
  pll_utree_t *tree = LibpllParsers::readNewickFromStr(geneTreeStr);
  int nodes = tree->tip_count + tree->inner_count;
  for (int i = 0; i < nodes; ++i) {
    auto node = tree->nodes[i];
    if (node->next) 
      continue;
    string label(node->label);
    string species;
    string gene;
    auto pos = label.find_first_of('_');
    species = label.substr(0, pos);
    gene = label.substr(pos + 1);
    Logger::info << gene << " " << species << endl;
    _map[gene] = species;
  }   
}

