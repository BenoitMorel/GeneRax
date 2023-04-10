#include "Families.hpp"
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>
#include <IO/FileSystem.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/LibpllException.hpp>
#include <algorithm>
#include <trees/PLLRootedTree.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <maths/Random.hpp>

enum FamilyErrorCode {
  ERROR_OK = 0,
  ERROR_READ_ALIGNMENT,
  ERROR_INVALID_LABEL,
  ERROR_READ_GENE_TREE,
  ERROR_GENE_TREE_POLYTOMY,
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
  case ERROR_GENE_TREE_POLYTOMY: return "Gene tree is not strictly binary";
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
  bool geneTreeProvided = family.startingGeneTree != "__random__" && family.startingGeneTree.size();
  bool mappingFileProvided = family.mappingFile.size();
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
    } catch (const LibpllException &e) {
      std::cerr << e.what() << std::endl;
    } catch(...) {
      return  ERROR_READ_ALIGNMENT;
    }
    if (alignmentLabels.size() < 3) {
      return ERROR_NOT_ENOUGH_GENES;
    }
  }
  if (mappingFileProvided) {
    if (!FileSystem::exists(family.mappingFile)) {
      return ERROR_MAPPING_FILE_EXISTENCE;
    }
  }
  // gene tree. Only check if one is given!
  if (geneTreeProvided) {
    if (!FileSystem::exists(family.startingGeneTree)) {
      return ERROR_GENE_TREE_FILE_EXISTENCE;
    }  
    std::unique_ptr<PLLUnrootedTree> geneTree; 
    try {
      geneTree = std::make_unique<PLLUnrootedTree>(
        family.startingGeneTree);
    } catch (const LibpllException &e) {
      std::cerr << e.what() << std::endl;
    } catch (...) {}
    if (!geneTree) {
      return ERROR_READ_GENE_TREE;
    } else if (!geneTree->isBinary()) {
      return ERROR_GENE_TREE_POLYTOMY;
    } else {
      geneTreeLabels = geneTree->getLabels();
      if (geneTreeLabels.size() < 3) {
        return ERROR_NOT_ENOUGH_GENES;
      }
      if (alignmentLabels.size() && alignmentLabels != geneTreeLabels) {
        std::cerr << "Absent in the gene tree: ";
        for (auto &l: alignmentLabels) {
          if (geneTreeLabels.find(l) == geneTreeLabels.end()) {
            std::cerr << l << " ";
          }
        }
        std::cerr << std::endl;
        std::cerr << "Absent in the alignment: ";
        for (auto &l: geneTreeLabels) {
          if (alignmentLabels.find(l) == alignmentLabels.end()) {
            std::cerr << l << " ";
          }
        }
        std::cerr << std::endl;
        std::cerr << std::endl;
        
        return ERROR_GENE_TREE_SEQUENCES_MISMATCH; 
      }
    }
    if (alignmentLabels.size() == 0 && geneTreeLabels.size() == 0) {
      return ERROR_NO_ALI_NO_TREE;
    }
  }
  GeneSpeciesMapping mapping;
  if (mappingFileProvided || geneTreeProvided) {
    mapping.fill(family.mappingFile, family.startingGeneTree);
  } else {
    mapping.fillFromGeneLabels(alignmentLabels);
  }
  if (speciesTreeLabels.size()) {
    auto &labels = (alignmentLabels.size() ? alignmentLabels : geneTreeLabels);
    if (!mapping.check(labels, speciesTreeLabels)) {
      return ERROR_MAPPING_MISMATCH;
    }
  }
  return ERROR_OK;
}

void Family::filterFamilies(Families &families, const std::string &speciesTreeFile, bool checkAlignments, bool checkSpeciesTree)
{
  ParallelContext::barrier();
  // at the end of this function, different ranks will have
  // a different rand state, so we save a seed
  auto consistentSeed = Random::getInt(); 
  Families copy = families;
  unsigned int initialFamilySize = static_cast<unsigned int>(copy.size());
  families.clear();
  std::unordered_set<std::string> speciesTreeLabels;
  unsigned int invalid = 0;
  corax_rtree_t *speciesTree = 0;
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
    corax_rtree_destroy(speciesTree, 0);
  }
  Random::setSeed(consistentSeed);
}



  
void Family::printStats(Families &families, 
    const std::string &speciesTreeFile, 
    const std::string &coverageFile,
    const std::string &fractionMissingFile)
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

  unsigned int maxGenes = 0;
  for (auto &species: speciesLabels) {
    maxGenes = std::max(maxGenes, 
        perSpeciesGenes[species]);
  }
  ParallelOfstream osFM(fractionMissingFile, true);
  for (auto &species: speciesLabels) {
    unsigned int genes = perSpeciesGenes[species];
    double fm = double(maxGenes - genes) / double(maxGenes);
    osFM << species << " " << fm << std::endl;
  }


  Logger::timed << "Input data information:" << std::endl;
  Logger::info << "- Number of gene families: " << families.size() << std::endl;
  Logger::info << "- Number of species: " << speciesTree.getLeafNumber() << std::endl;
  Logger::info << "- Total number of genes: " << totalGeneNumber << std::endl;
  Logger::info << "- Average number of genes per family: " << totalGeneNumber / families.size() << std::endl;
  Logger::info << "- Maximum number of genes per family: " << maxGeneNumber << std::endl;
  Logger::info << "- Species covered with the smallest family coverage: \"" << minCoveredSpecies << "\" (covered by " << minSpeciesCoverage << "/" << families.size() << " families)" << std::endl;
  Logger::info << "- Average (over species) species family coverage: " << totalSpeciesCoverage / speciesLabels.size() << std::endl;
  Logger::info << std::endl;
}


