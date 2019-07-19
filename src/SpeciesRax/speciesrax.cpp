#include "SpeciesRaxArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <limits>
#include <IO/FileSystem.hpp>
#include <sstream>
#include <optimizers/SpeciesTreeOptimizer.hpp>
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
  
  SpeciesTreeOptimizer speciesTreeOptimizer(arguments.speciesTree, initialFamilies, recModel, arguments.output);
  speciesTreeOptimizer.sprSearch(1);
  speciesTreeOptimizer.rootExhaustiveSearch();
  speciesTreeOptimizer.ratesOptimization();
  speciesTreeOptimizer.sprSearch(2);
  speciesTreeOptimizer.rootExhaustiveSearch();
  speciesTreeOptimizer.ratesOptimization();
  speciesTreeOptimizer.sprSearch(3);

  Logger::timed << "End of the run" << std::endl;
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

