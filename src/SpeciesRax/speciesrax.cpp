#include "SpeciesRaxArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <limits>
#include <IO/FileSystem.hpp>
#include <sstream>
#include <trees/SpeciesTree.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <memory>

void initFolders(const std::string &output, Families &families) 
{
  std::string results = FileSystem::joinPaths(output, "results");
  FileSystem::mkdir(results, true);
  for (auto &family: families) {
    FileSystem::mkdir(FileSystem::joinPaths(results, family.name), true);
  }
}


int internal_main(int argc, char** argv, void* comm)
{
  // the order is very important
  ParallelContext::init(comm); 
  Logger::init();
  SpeciesRaxArguments arguments(argc, argv);
  FileSystem::mkdir(arguments.output, true);
  Logger::initFileOutput(FileSystem::joinPaths(arguments.output, "speciesrax"));
  srand(arguments.seed);
  arguments.printCommand();
  arguments.printSummary();
  
  RecModel recModel = arguments.reconciliationModel;

  Families initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << std::endl;
  initFolders(arguments.output, initialFamilies);
  
  
  PerCoreGeneTrees geneTrees(initialFamilies); 
  std::unique_ptr<SpeciesTree> speciesTree;
  Logger::info << "species tree: " << arguments.speciesTree << std::endl;
  if (arguments.speciesTree == "random") {
    speciesTree = std::make_unique<SpeciesTree>(initialFamilies);
    DTLRates rates(0.1, 0.2, 0.1);
    speciesTree->setRates(rates);
  } else {
    speciesTree = std::make_unique<SpeciesTree>(arguments.speciesTree);
    SpeciesTreeOptimizer::ratesOptimization(*speciesTree, geneTrees, recModel);
  }
  speciesTree->saveToFile(FileSystem::joinPaths(arguments.output, "starting_species_tree.newick"));
  Logger::info << *speciesTree << std::endl;
  SpeciesTreeOptimizer::sprSearch(*speciesTree, geneTrees, recModel, 1);
  SpeciesTreeOptimizer::rootExhaustiveSearch(*speciesTree, geneTrees, recModel);
  SpeciesTreeOptimizer::ratesOptimization(*speciesTree, geneTrees, recModel);
  SpeciesTreeOptimizer::sprSearch(*speciesTree, geneTrees, recModel, 2);
  SpeciesTreeOptimizer::rootExhaustiveSearch(*speciesTree, geneTrees, recModel);
  SpeciesTreeOptimizer::ratesOptimization(*speciesTree, geneTrees, recModel);
  SpeciesTreeOptimizer::sprSearch(*speciesTree, geneTrees, recModel, 2);
  SpeciesTreeOptimizer::ratesOptimization(*speciesTree, geneTrees, recModel);
  SpeciesTreeOptimizer::sprSearch(*speciesTree, geneTrees, recModel, 3);

  Logger::timed << "End of the run" << std::endl;
  speciesTree->saveToFile(FileSystem::joinPaths(arguments.output, "inferred_species_tree.newick"));
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

