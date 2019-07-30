#include "SpeciesRaxArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <routines/GeneRaxSlaves.hpp>
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
  std::string proposals = FileSystem::joinPaths(output, "proposals");
  FileSystem::mkdir(results, true);
  FileSystem::mkdir(proposals, true);
  for (auto &family: families) {
    FileSystem::mkdir(FileSystem::joinPaths(results, family.name), true);
    FileSystem::mkdir(FileSystem::joinPaths(proposals, family.name), true);
  }
}


int speciesrax_main(int argc, char** argv, void* comm)
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
  
  SpeciesTreeOptimizer speciesTreeOptimizer(arguments.speciesTree, initialFamilies, recModel, arguments.output, argv[0]);
  for (unsigned int radius = 1; radius <= arguments.fastRadius; ++radius) {
    speciesTreeOptimizer.ratesOptimization();
    speciesTreeOptimizer.sprSearch(radius, false);
    speciesTreeOptimizer.rootExhaustiveSearch(false);
  }
  for (unsigned int radius = 1; radius <= arguments.slowRadius; ++radius) {
    speciesTreeOptimizer.advancedRatesOptimization(1);
    speciesTreeOptimizer.sprSearch(radius, true);
    speciesTreeOptimizer.rootExhaustiveSearch(true);
  } 
  speciesTreeOptimizer.rootExhaustiveSearch(true);
  Logger::timed << "End of the run" << std::endl;
  ParallelContext::finalize();
  return 0;
}



/**
 * GeneRax can call itself with the scheduler. This function
 * decides whether we are in the main called by the user (generax_main)
 * or if we are called by the scheduler to execute some 
 * intermediate step (static_scheduled_main)
 */
int internal_main(int argc, char** argv, void* comm)
{
  if (GeneRaxSlaves::is_slave(argc, argv)) {
    int slaveComm = -1; 
    return static_scheduled_main(argc, argv, &slaveComm);
  } else {
    return speciesrax_main(argc, argv, comm);
  }
}

int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}


