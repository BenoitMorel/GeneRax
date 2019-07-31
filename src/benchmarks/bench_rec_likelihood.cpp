#include <trees/SpeciesTree.hpp>
#include <cassert>
#include <memory>
#include <families/Families.hpp>
#include <IO/FileSystem.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <optimizers/SpeciesTreeOptimizer.hpp>
#include <routines/GeneRaxSlaves.hpp>

void generateFakeAlignment(const std::vector<std::string> &taxa, const std::string &outputFile)
{
  std::ofstream os(outputFile);
  for (auto &taxon: taxa) { 
    os << ">" << taxon << std::endl << "ACGT" << std::endl;
  }
}

void generateRandomFamilies(const std::string &outputDirectory,
    std::unordered_set<std::string> &speciesLabels,
    unsigned int familiesNumber,
    unsigned int geneTaxaNumber,
    Families &families)
{
  std::vector<std::string> vectorSpeciesLabels;
  for (auto &str: speciesLabels) {
    vectorSpeciesLabels.push_back(str);
  }
  for (unsigned int i = 0; i < familiesNumber; ++i) {
    FamilyInfo family;
    family.name = std::string("family_" + std::to_string(i));
    family.startingGeneTree = FileSystem::joinPaths(outputDirectory, family.name + ".newick");
    family.alignmentFile = FileSystem::joinPaths(outputDirectory, family.name + ".fasta");
    std::vector<std::string> geneTaxa;
    for (unsigned int j = 0; j < geneTaxaNumber; ++j) {
      std::string taxa = vectorSpeciesLabels[rand() % vectorSpeciesLabels.size()];
      taxa += "_gene" + std::to_string(j);
      geneTaxa.push_back(taxa);
    }
    generateFakeAlignment(geneTaxa, family.alignmentFile);
    LibpllEvaluation::createAndSaveRandomTree(family.alignmentFile, family.libpllModel, family.startingGeneTree);
    families.push_back(family);
  }
}

void bench(const std::string &speciesTree, 
    Families &families, 
    RecModel recModel,
    unsigned int fastRadius,
    const std::string &outputDirectory) 
{
  SpeciesTreeOptimizer speciesTreeOptimizer(speciesTree, families, recModel, outputDirectory, "");
  for (unsigned int radius = 1; radius <= fastRadius; ++radius) {
    speciesTreeOptimizer.ratesOptimization();
    speciesTreeOptimizer.sprSearch(radius, false);
    speciesTreeOptimizer.rootExhaustiveSearch(false);

  }
}

// hack to compile the project
void plop()
{
  int slaveComm;
  static_scheduled_main(0, 0, &slaveComm);
}

int main(int, char**)
{
  // Input parameters
  unsigned int speciesNumber = 15;
  unsigned int familiesNumber = 10;
  unsigned int geneTaxaNumber = 20;
  unsigned int fastRadius = 3;
  RecModel recModel = UndatedDTL;
  
  
  
  srand(42);
  std::string outputDirectory("tmp");
  std::string dataDirectory = FileSystem::joinPaths(outputDirectory, "data");
  std::string runDirectory = FileSystem::joinPaths(outputDirectory, "run");

  FileSystem::mkdir(outputDirectory, true);
  FileSystem::mkdir(dataDirectory, true);
  FileSystem::mkdir(runDirectory, true);
  std::unordered_set<std::string> speciesLabels;
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    speciesLabels.insert(std::string("S" + std::to_string(i)));
  }
  SpeciesTree speciesTree(speciesLabels);
  std::string speciesTreeFile = FileSystem::joinPaths(dataDirectory, "speciesTree.newick");
  speciesTree.saveToFile(speciesTreeFile);
  Families families;
  generateRandomFamilies(dataDirectory, 
      speciesLabels,
      familiesNumber,
      geneTaxaNumber,
      families);      
  bench(speciesTreeFile, families, recModel, fastRadius, runDirectory);

  return 0;
}

