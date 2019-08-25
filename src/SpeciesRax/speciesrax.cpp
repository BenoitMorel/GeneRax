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

void initFamilies(const std::string &output, Families &families) 
{
  std::string results = FileSystem::joinPaths(output, "results");
  std::string proposals = FileSystem::joinPaths(output, "proposals");
  FileSystem::mkdir(results, true);
  FileSystem::mkdir(proposals, true);
  for (auto &family: families) {
    std::string familyResultsDir = FileSystem::joinPaths(results, family.name);
    FileSystem::mkdir(familyResultsDir, true);
    FileSystem::mkdir(FileSystem::joinPaths(proposals, family.name), true);
    std::string inputGeneTree = family.startingGeneTree;
    std::string outputGeneTree = FileSystem::joinPaths(familyResultsDir, family.name + ".newick");
    FileSystem::copy(inputGeneTree, outputGeneTree, true);
    family.startingGeneTree = outputGeneTree;
  }
  ParallelContext::barrier();
}


void simpleSearch(SpeciesRaxArguments &arguments, char ** argv)
{
  RecModel recModel = arguments.reconciliationModel;
  Families initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << std::endl;
  initFamilies(arguments.output, initialFamilies);
  SpeciesTreeOptimizer speciesTreeOptimizer(arguments.speciesTree, initialFamilies, UndatedDL, arguments.output, argv[0]);
  speciesTreeOptimizer.setPerSpeciesRatesOptimization(arguments.perSpeciesDTLRates); 
  for (unsigned int radius = 1; radius <= arguments.fastRadius; ++radius) {
    if (radius == arguments.fastRadius) {
      speciesTreeOptimizer.setModel(recModel);
    }
    speciesTreeOptimizer.ratesOptimization();
    speciesTreeOptimizer.sprSearch(radius, false);
    speciesTreeOptimizer.rootExhaustiveSearch(false);
    Logger::info << "RecLL = " << speciesTreeOptimizer.getReconciliationLikelihood() << std::endl;
  }
  if (arguments.slowRadius) {
    Logger::info << "FIRST SLOW RADIUS SEARCH, WITHOUT TRANSFER " << std::endl;
    speciesTreeOptimizer.setModel(UndatedDL); 
    speciesTreeOptimizer.sprSearch(arguments.fastRadius, true);
    speciesTreeOptimizer.rootExhaustiveSearch(true);
    if (recModel != UndatedDL) {
      Logger::info << "SECOND SLOW RADIUS SEARCH WITH TRANSFERS " << std::endl;
      speciesTreeOptimizer.setModel(recModel); 
      speciesTreeOptimizer.sprSearch(arguments.fastRadius, true);
    }
    Logger::info << "Joint LL = " << speciesTreeOptimizer.computeLikelihood(true) << std::endl;
  }
  if (arguments.finalGeneRadius) {
    Logger::info << "Final gene tree optimization step, with radius " << arguments.finalGeneRadius << std::endl;
    speciesTreeOptimizer.optimizeGeneTrees(arguments.finalGeneRadius, true);
    Logger::info << "Joint LL = " << speciesTreeOptimizer.computeLikelihood(true) << std::endl;
  }
  speciesTreeOptimizer.saveCurrentSpeciesTree();
}

void subsampleSearch(SpeciesRaxArguments &arguments, char ** argv)
{
  RecModel recModel = arguments.reconciliationModel;
  Families initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << std::endl;
  initFamilies(arguments.output, initialFamilies);
  SpeciesTreeOptimizer speciesTreeOptimizer(arguments.speciesTree, initialFamilies, recModel, arguments.output, argv[0]);
  speciesTreeOptimizer.setPerSpeciesRatesOptimization(arguments.perSpeciesDTLRates); 
  unsigned int sampleSize = initialFamilies.size() / 3;
  unsigned int sampleNumber = 5;
  unsigned int iterations = 1;
  for (unsigned int it = 0; it < iterations; ++it) {
    std::unordered_set<std::string> speciesIds;
    for (unsigned int i = 0; i < sampleNumber; ++i) {
      std::string speciesId = std::string("sample_") + std::to_string(it) + std::string("_") + std::to_string(i);
      speciesTreeOptimizer.inferSpeciesTreeFromSamples(sampleSize, speciesId);
      speciesIds.insert(speciesId);
    }
    for (unsigned int i =0; i < 1; ++i) {
      speciesTreeOptimizer.optimizeGeneTreesFromSamples(speciesIds, std::string("opt_sub_genetrees_") + std::to_string(it) + std::string("_") + std::to_string(i));
    }
    Logger::info << "RecLL = " << speciesTreeOptimizer.getReconciliationLikelihood() << std::endl;
  }
  
  
  for (unsigned int radius = 1; radius <= arguments.fastRadius; ++radius) {
    if (radius == arguments.fastRadius) {
      speciesTreeOptimizer.setModel(recModel);
    }
    speciesTreeOptimizer.ratesOptimization();
    speciesTreeOptimizer.sprSearch(radius, false);
    speciesTreeOptimizer.rootExhaustiveSearch(false);
    Logger::info << "RecLL = " << speciesTreeOptimizer.getReconciliationLikelihood() << std::endl;
    Logger::info << "Joint LL = " << speciesTreeOptimizer.computeLikelihood(true) << std::endl;
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
  
  switch (arguments.strategy) {
  case SIMPLE_SEARCH:
    simpleSearch(arguments, argv);
    break;
  case SUBSAMPLE_SEARCH:
    subsampleSearch(arguments, argv);
    break;
  }
  Logger::timed << "Output in " << arguments.output << std::endl;
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


