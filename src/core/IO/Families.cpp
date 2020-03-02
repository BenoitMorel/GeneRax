#include "Families.hpp"
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <IO/FileSystem.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <algorithm>
#include <trees/PLLRootedTree.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>

enum FamilyErrorCode {
  ERROR_OK = 0,
  ERROR_READ_ALIGNMENT,
  ERROR_INVALID_LABEL,
  ERROR_READ_GENE_TREE,
  ERROR_NOT_ENOUGH_GENES,
  ERROR_GENE_TREE_SEQUENCES_MISMATCH,
  ERROR_ALIGNEMENT_FILE_EXISTENCE,
  ERROR_GENE_TREE_FILE_EXISTENCE,
  ERROR_MAPPING_FILE_EXISTENCE,
  ERROR_MAPPING_MISMATCH,
  ERROR_NO_ALI_NO_TREE
};

static std::string getErrorMessage(FamilyErrorCode error) {
  switch(error) {
  case ERROR_OK: 
    assert(0); 
    return ""; 
  case ERROR_READ_ALIGNMENT: return "Cannot read alignment file (file exists but is invalid)";
  case ERROR_INVALID_LABEL: return "Some labels in the MSAs contain invalid characters";
  case ERROR_READ_GENE_TREE: return "Cannot read starting gene tree file (file exists but is invalid)";
  case ERROR_NOT_ENOUGH_GENES : return "The gene family should contain at least 3 sequences";
  case ERROR_GENE_TREE_SEQUENCES_MISMATCH : return "Mismatch between gene tree leaves and sequence labels";
  case ERROR_ALIGNEMENT_FILE_EXISTENCE : return "Input alignment file does not exist";
  case ERROR_GENE_TREE_FILE_EXISTENCE : return "Input starting gene tree file does not exist";
  case ERROR_MAPPING_FILE_EXISTENCE : return "A mapping file was given but does not exist";
  case ERROR_MAPPING_MISMATCH : return "Failed to map genes and species";
  case ERROR_NO_ALI_NO_TREE : return "This family does not have any starting gene tree nor alignment";
  }
  return "";
}

static FamilyErrorCode filterFamily(const FamilyInfo &family, const std::unordered_set<std::string> &speciesTreeLabels, bool checkAlignments)
{
  std::unordered_set<std::string> alignmentLabels;
  std::unordered_set<std::string> geneTreeLabels;
  // sequences
  if (checkAlignments) {
    try {
      if (!FileSystem::exists(family.alignmentFile)) {
        return ERROR_ALIGNEMENT_FILE_EXISTENCE;
      }
      if (!LibpllParsers::fillLabelsFromAlignment(family.alignmentFile, family.libpllModel, alignmentLabels)) {
        return ERROR_READ_ALIGNMENT;
      }
      if (!LibpllParsers::areLabelsValid(alignmentLabels)) {
        return ERROR_INVALID_LABEL;
      }
    } catch (LibpllException e) {
      std::cerr << e.what() << std::endl;
    } catch(...) {
      return  ERROR_READ_ALIGNMENT;
    }
    if (alignmentLabels.size() < 3) {
      return ERROR_NOT_ENOUGH_GENES;
    }
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
      if (geneTreeLabels.size() < 3) {
        return ERROR_NOT_ENOUGH_GENES;
      }
      if (alignmentLabels.size() && alignmentLabels != geneTreeLabels) {
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
    if (alignmentLabels.size() == 0 && geneTreeLabels.size() == 0) {
      return ERROR_NO_ALI_NO_TREE;
    }
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    if (speciesTreeLabels.size()) {
      auto &labels = (alignmentLabels.size() ? alignmentLabels : geneTreeLabels);
      if (!mapping.check(labels, speciesTreeLabels)) {
        return ERROR_MAPPING_MISMATCH;
      }
    }
  }
  return ERROR_OK;
}

void Family::filterFamilies(Families &families, const std::string &speciesTreeFile, bool checkAlignments, bool checkSpeciesTree)
{
  ParallelContext::barrier();
  // at the end of this function, different ranks will have
  // a different rand state, so we save a seed
  auto consistentSeed = rand(); 
  Families copy = families;
  unsigned int initialFamilySize = static_cast<unsigned int>(copy.size());
  families.clear();
  std::unordered_set<std::string> speciesTreeLabels;
  unsigned int invalid = 0;
  pll_rtree_t *speciesTree = 0;
  if (checkSpeciesTree) {
    if (!FileSystem::exists(speciesTreeFile)) {
      Logger::info << "[Error] Species tree file does not exist (" << speciesTreeFile << ")" << std::endl;
      ParallelContext::abort(10);
    }
    try {
      speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile);
    } catch (...) {}
    if (!speciesTree) {
      Logger::info << "[Error] Cannot parse species tree file (" << speciesTreeFile << ")" << std::endl;
      ParallelContext::abort(10);
    }
    LibpllParsers::fillLeavesFromRtree(speciesTree, speciesTreeLabels);
  } 
  std::vector<unsigned int> localErrors((initialFamilySize - 1 ) / ParallelContext::getSize() + 1, 99999999);
  std::vector<unsigned int> errors;
  for (auto i = ParallelContext::getBegin(initialFamilySize); i < ParallelContext::getEnd(initialFamilySize); i ++) {
    auto &family = copy[i];
    localErrors[i - ParallelContext::getBegin(initialFamilySize)]  = filterFamily(family, speciesTreeLabels, checkAlignments); 
  }
  ParallelContext::concatenateUIntVectors(localErrors, errors);
  errors.erase(remove(errors.begin(), errors.end(), 99999999), errors.end());
  Logger::info << std::endl;
  ParallelContext::barrier();
  for (unsigned int i = 0; i < initialFamilySize; ++i) {
    auto error = errors[i];
    auto &family = copy[i];
    if (ERROR_OK == error) {
      families.push_back(family);
    } else {
      Logger::info << "Error in family " << family.name 
        << ": "  << getErrorMessage(static_cast<FamilyErrorCode>(error)) << std::endl;
      invalid++;
    }
  }
  if (invalid) {
    Logger::info << "WARNING!!! Found " << invalid 
      << " invalid families (they will be discarded from the analysis)" << std::endl;
  }  
  if (speciesTree) {
    pll_rtree_destroy(speciesTree, 0);
  }
  srand(consistentSeed);
}



  
void Family::printStats(Families &families, const std::string &speciesTreeFile, const std::string &coverageFile)
{
 
  PLLRootedTree speciesTree(speciesTreeFile);

  std::unordered_map<std::string, unsigned int> perSpeciesGenes;
  std::unordered_map<std::string, unsigned int> perSpeciesCoveringFamilies;
  unsigned int totalGeneNumber = 0;
  unsigned int maxGeneNumber = 0;
  unsigned int totalSpeciesCoverage = 0;
  unsigned int minSpeciesCoverage = 999999999;
  std::string minCoveredSpecies;

  auto speciesLabels = speciesTree.getLabels(true);

  for (const auto &species: speciesLabels) {
    perSpeciesGenes.insert({species, 0});
    perSpeciesCoveringFamilies.insert({species, 0});
  }

  PerCoreGeneTrees geneTrees(families);
  for (auto &tree: geneTrees.getTrees()) {
    for (auto species: tree.mapping.getCoveredSpecies()) {
      perSpeciesCoveringFamilies[species]++;
    }
    for (auto geneSpecies: tree.mapping.getMap()) {
      perSpeciesGenes[geneSpecies.second]++;
    }
    unsigned int geneNumber = tree.mapping.getMap().size();
    totalGeneNumber += geneNumber;
    maxGeneNumber = std::max(geneNumber, maxGeneNumber);
  }


  // gather parallel values
  ParallelContext::sumUInt(totalGeneNumber);
  ParallelContext::maxUInt(maxGeneNumber);
  for (const auto &species: speciesLabels) {
    ParallelContext::sumUInt(perSpeciesCoveringFamilies[species]);
    ParallelContext::sumUInt(perSpeciesGenes[species]);
    totalSpeciesCoverage += perSpeciesCoveringFamilies[species];
    if (minSpeciesCoverage > perSpeciesCoveringFamilies[species]) {
      minSpeciesCoverage = perSpeciesCoveringFamilies[species];
      minCoveredSpecies = species;
    }
  }

  ParallelOfstream os(coverageFile, true);
  typedef std::pair<double, std::string> Coverage;
  std::vector<Coverage> coverageVector;
  for (const auto &species: speciesLabels) {
    double d = static_cast<double>(perSpeciesCoveringFamilies[species]);
    d /= static_cast<double>(families.size());
    coverageVector.push_back({d, species});
  }
  std::sort(coverageVector.begin(), coverageVector.end());
  os << "SPECIES: FAMILY_COVERAGE" << std::endl;
  for (auto &coverage: coverageVector) {
    os << coverage.second << ": " << coverage.first << std::endl;
  }

  Logger::timed << "Input data information:" << std::endl;
  Logger::info << "- Number of gene families: " << families.size() << std::endl;
  Logger::info << "- Number of species: " << speciesTree.getLeavesNumber() << std::endl;
  Logger::info << "- Total number of genes: " << totalGeneNumber << std::endl;
  Logger::info << "- Average number of genes per family: " << totalGeneNumber / families.size() << std::endl;
  Logger::info << "- Maximum number of genes per family: " << maxGeneNumber << std::endl;
  Logger::info << "- Species covered with the smallest family coverage: \"" << minCoveredSpecies << "\" (covered by " << minSpeciesCoverage << "/" << families.size() << " families)" << std::endl;
  Logger::info << "- Average (over species) species family coverage: " << totalSpeciesCoverage / speciesLabels.size() << std::endl;
  Logger::info << std::endl;
}


