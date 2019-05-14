#include "GeneSpeciesMapping.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>

void GeneSpeciesMapping::fill(const std::string &mappingFile, const std::string &geneTreeStr) 
{
  if (mappingFile.size()) {
    buildFromMappingFile(mappingFile);
  } else {
    buildFromTrees(geneTreeStr);
  }
}

void GeneSpeciesMapping::buildFromMappingFile(const std::string &mappingFile)
{
  std::ifstream f(mappingFile);
  std::string line;
  getline(f, line);
  f.seekg(0, ios::beg);
  if (line.find(":") != std::string::npos) {
    buildFromPhyldogMapping(f);
  } else {
    buildFromTreerecsMapping(f);
  }
}
  

void GeneSpeciesMapping::buildFromPhyldogMapping(std::ifstream &f)
{
  /*
    species1:gene1;gene2;gene3
    species2:gene4;gene5
  */
  std::string line;
  while (getline(f, line)) {
    stringstream ss(line);
    std::string species;
    std::string gene;
    getline(ss, species, ':');
    while(getline(ss, gene, ';')) {
      _map[gene] = species;
    }
  }
}

void GeneSpeciesMapping::buildFromTreerecsMapping(std::ifstream &f)
{
  /*
    gene1 species1
    gene2 species2
    gene3 species1
  */
  std::string line;
  while (getline(f, line)) {
    stringstream ss(line);
    std::string species;
    std::string gene;
    ss >> gene;
    ss >> species;
    _map[gene] = species;
  }
}

void GeneSpeciesMapping::buildFromTrees(const std::string &geneTreeStr)
{
  pll_utree_t *tree = LibpllParsers::readNewickFromStr(geneTreeStr);
  int nodes = tree->tip_count + tree->inner_count;
  for (int i = 0; i < nodes; ++i) {
    auto node = tree->nodes[i];
    if (node->next) 
      continue;
    std::string label(node->label);
    std::string species;
    std::string gene;
    auto pos = label.find_first_of('_');
    species = label.substr(0, pos);
    gene = label.substr(pos + 1);
    Logger::info << gene << " " << species << std::endl;
    _map[gene] = species;
  }   
}

