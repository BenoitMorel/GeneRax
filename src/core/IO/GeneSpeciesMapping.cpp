#include "GeneSpeciesMapping.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>
extern "C" {
#include <pll.h>
}

void GeneSpeciesMapping::fill(const std::string &mappingFile, const std::string &geneTreeStr) 
{
  if (mappingFile.size()) {
    buildFromMappingFile(mappingFile);
  } else {
    buildFromTrees(geneTreeStr);
  }
}
  
void GeneSpeciesMapping::fill(const GeneSpeciesMapping &mapping)
{
  auto &m = mapping.getMap();
  for (auto &pair: m) {
    _map[pair.first] = pair.first;
  }
}
  
bool GeneSpeciesMapping::check(pll_utree_t *geneTree, pll_rtree_t *speciesTree)
{
  std::unordered_set<std::string> geneLeaves;
  std::unordered_set<std::string> speciesLeaves;
  LibpllParsers::fillLeavesFromUtree(geneTree, geneLeaves); 
  LibpllParsers::fillLeavesFromRtree(speciesTree, speciesLeaves); 
  bool ok = true;
  for (auto pair: getMap()) {
    auto &gene = pair.first;
    auto &species = pair.second;
    if (geneLeaves.find(gene) == geneLeaves.end()) {
      Logger::error << "[Error] Invalid mapping " << gene << "<->" << species << ": can't find the gene " << gene << " in the gene tree" << std::endl;
      ok = false;
    }
    if (speciesLeaves.find(species) == speciesLeaves.end()) {
      Logger::error << "[Error] Invalid mapping " << gene << "<->" << species << ": can't find the species " << species << " in the species tree" << std::endl;
      ok = false;
    }
  }
  for (auto &gene: geneLeaves) {
    auto speciesIt = getMap().find(gene);
    if (speciesIt == getMap().end()) {
      Logger::error << "[Error] Gene tree leaf " << gene << " is not mapped to any species" << std::endl;
      ok = false;
    }
  }
  return ok;
}


void GeneSpeciesMapping::buildFromMappingFile(const std::string &mappingFile)
{
  std::ifstream f(mappingFile);
  std::string line;
  getline(f, line);
  f.seekg(0, std::ios::beg);
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
    std::stringstream ss(line);
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
    std::stringstream ss(line);
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
  auto nodes = tree->tip_count + tree->inner_count;
  for (unsigned int i = 0; i < nodes; ++i) {
    auto node = tree->nodes[i];
    if (node->next) 
      continue;
    std::string label(node->label);
    std::string species;
    std::string gene;
    auto pos = label.find_first_of('_');
    species = label.substr(0, pos);
    gene = label; //label.substr(pos + 1);
    _map[gene] = species;
  }   
}

