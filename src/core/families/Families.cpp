#include "Families.hpp"
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <algorithm>

enum FamilyErrorCode {
  ERROR_OK = 0,
  ERROR_READ_ALIGNMENT,
  ERROR_READ_GENE_TREE,
  ERROR_NOT_ENOUGH_GENES,
  ERROR_GENE_TREE_SEQUENCES_MISMATCH,
  ERROR_ALIGNEMENT_FILE_EXISTENCE,
  ERROR_GENE_TREE_FILE_EXISTENCE,
  ERROR_MAPPING_FILE_EXISTENCE,
  ERROR_MAPPING_MISMATCH
};

std::string getErrorMessage(FamilyErrorCode error) {
  switch(error) {
  case ERROR_OK: assert(0); return ""; 
  case ERROR_READ_ALIGNMENT: return "Cannot read alignment file (file exists but is invalid)";
  case ERROR_READ_GENE_TREE: return "Cannot read starting gene tree file (file exists but is invalid)";
  case ERROR_NOT_ENOUGH_GENES : return "The gene family should contain at least 3 sequences";
  case ERROR_GENE_TREE_SEQUENCES_MISMATCH : return "Mismatch between gene tree leaves and sequence labels";
  case ERROR_ALIGNEMENT_FILE_EXISTENCE : return "Input alignment file does not exist";
  case ERROR_GENE_TREE_FILE_EXISTENCE : return "Input starting gene tree file does not exist";
  case ERROR_MAPPING_FILE_EXISTENCE : return "A mapping file was given but does not exist";
  case ERROR_MAPPING_MISMATCH : return "Failed to map genes and species";
  }
  return "";
}

static FamilyErrorCode filterFamily(const FamilyInfo &family, const std::unordered_set<std::string> &speciesTreeLabels)
{
  std::unordered_set<std::string> alignmentLabels;
  // sequences
  try {
    if (!FileSystem::exists(family.alignmentFile)) {
      return ERROR_ALIGNEMENT_FILE_EXISTENCE;
    }
    if (!LibpllParsers::fillLabelsFromAlignment(family.alignmentFile, family.libpllModel, alignmentLabels)) {
      return ERROR_READ_ALIGNMENT;
    }
  } catch (LibpllException e) {
    std::cerr << e.what() << std::endl;
  } catch(...) {
    return  ERROR_READ_ALIGNMENT;
  }
  if (alignmentLabels.size() < 3) {
    return ERROR_NOT_ENOUGH_GENES;
  }
  if (family.mappingFile.size()) {
    if (!FileSystem::exists(family.mappingFile)) {
      return ERROR_MAPPING_FILE_EXISTENCE;
    }
  }
  // gene tree. Only check if one is given!
  if (family.startingGeneTree != "__random__" && family.startingGeneTree.size()) {
    if (!FileSystem::exists(family.startingGeneTree)) {
      return ERROR_GENE_TREE_FILE_EXISTENCE;
    }  
    std::unordered_set<std::string> geneTreeLabels;
    pll_utree_t * utree = 0;
    try {
      utree = LibpllParsers::readNewickFromFile(family.startingGeneTree);
    } catch (LibpllException e) {
      std::cerr << e.what() << std::endl;
    } catch (...) {}
    if (!utree) {
      return ERROR_READ_GENE_TREE;
    } else {
      LibpllParsers::fillLeavesFromUtree(utree, geneTreeLabels);
      pll_utree_destroy(utree, 0);
      if (alignmentLabels != geneTreeLabels) {
        for (auto &l: alignmentLabels) {
          std::cerr << l << " ";
        }
        std::cerr << std::endl;
        for (auto &l: geneTreeLabels) {
          std::cerr << l << " ";
        }
        std::cerr << std::endl;
        std::cerr << std::endl;
        
        return ERROR_GENE_TREE_SEQUENCES_MISMATCH; 
      }
    }
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    if (!mapping.check(alignmentLabels, speciesTreeLabels)) {
      return ERROR_MAPPING_MISMATCH;
    }
  }
  return ERROR_OK;
}

void filterFamilies(Families &families, const std::string &speciesTreeFile)
{
  ParallelContext::barrier();
  Families copy = families;
  families.clear();
  std::unordered_set<std::string> speciesTreeLabels;
  unsigned int invalid = 0;
  if (!FileSystem::exists(speciesTreeFile)) {
    Logger::info << "[Error] Species tree file does not exist (" << speciesTreeFile << ")" << std::endl;
    ParallelContext::abort(10);
  }
  pll_rtree_t *speciesTree = 0;
  try {
    speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile);
  } catch (...) {}
  if (!speciesTree) {
    Logger::info << "[Error] Cannot parse species tree file (" << speciesTreeFile << ")" << std::endl;
    ParallelContext::abort(10);
  }
  LibpllParsers::fillLeavesFromRtree(speciesTree, speciesTreeLabels);
  
  std::vector<unsigned int> localErrors((copy.size() - 1 ) / ParallelContext::getSize() + 1, 99999999);
  std::vector<unsigned int> errors;
  for (auto i = ParallelContext::getBegin(copy.size()); i < ParallelContext::getEnd(copy.size()); i ++) {
    auto &family = copy[i];
    localErrors[i - ParallelContext::getBegin(copy.size())]  = filterFamily(family, speciesTreeLabels); 
  }
  ParallelContext::concatenateUIntVectors(localErrors, errors);
  errors.erase(remove(errors.begin(), errors.end(), 99999999), errors.end());
  Logger::info << std::endl;
  ParallelContext::barrier();
  for (unsigned int i = 0; i < copy.size(); ++i) {
    auto error = errors[i];
    auto &family = copy[i];
    if (ERROR_OK == error) {
      families.push_back(family);
    } else {
      Logger::info << "Error in family " << family.name << ": "  << getErrorMessage(static_cast<FamilyErrorCode>(error)) << std::endl;
      invalid++;
    }
  }
  if (invalid) {
    Logger::info << "WARNING!!! Found " << invalid << " invalid families (they will be discarded from the analysis)" << std::endl;
  }   
}

void duplicatesFamilies(const Families &families, Families &duplicatedFamilies, unsigned int factor)
{
  duplicatedFamilies.clear();
  for (auto &family: families) {
    for (unsigned int i = 0; i < factor; ++i) {
      FamilyInfo child = family;
      child.name += "_" + std::to_string(i);
      duplicatedFamilies.push_back(child);
    }
  }
}

double getLL(const FamilyInfo &family)
{
    std::ifstream is(family.statsFile);
    double libpllLL = 0.0;
    double recLL = 0.0;
    is >> libpllLL;
    is >> recLL;
    return libpllLL + recLL;
}
void contractFamilies(const Families &duplicatedFamilies, Families &families)
{
  unsigned int factor = duplicatedFamilies.size() / families.size();
  for (unsigned int f = 0; f < families.size(); ++f) {
    double bestLL = -999999999;
    unsigned int bestIndex = 0;
    unsigned int bestChunk = 0;
    for (unsigned int i = 0; i < factor; ++i) {
      auto &duplicatedFamily = duplicatedFamilies[f * factor + i];
      auto ll = getLL(duplicatedFamily);
      if (ll > bestLL) {
        bestIndex = f * factor + i;
        bestLL = ll;
        bestChunk = i;
      }
    }
    auto name = families[f].name;
    families[f] = duplicatedFamilies[bestIndex];
    families[f].name = name;
    Logger::info << "color " << families[f].color  <<  " chunk: " << bestChunk << std::endl;
  }

}


void splitInitialFamilies(const Families &families, std::vector<Families> &splitFamilies, unsigned int splitsNumber)
{
  splitFamilies.clear();
  splitFamilies.resize(splitsNumber);
  for (unsigned int i = 0; i < families.size(); ++i) {
    splitFamilies[i % splitsNumber].push_back(families[i]);
    splitFamilies[i % splitsNumber].back().color = i % splitsNumber;
    assert(families[i].name ==  splitFamilies[i % splitsNumber][i / splitsNumber].name);

  }
}


void mergeSplitFamilies(const std::vector<Families> &splitFamilies, Families &families, unsigned int splitsNumber)
{
  for (unsigned int i = 0; i < families.size(); ++i) {
    assert(families[i].name ==  splitFamilies[i % splitsNumber][i / splitsNumber].name);
    families[i]  = splitFamilies[i % splitsNumber][i / splitsNumber];
  }
}
