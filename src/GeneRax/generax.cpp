#include "GeneRaxArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <limits>
#include <trees/PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/DTLRates.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <routines/GeneRaxSlaves.hpp>
#include <parallelization/Scheduler.hpp>
#include <routines/RaxmlMaster.hpp>
#include <routines/GeneTreeSearchMaster.hpp>
#include <routines/Routines.hpp>

bool useSplitImplem() {
  return ParallelContext::getSize() > 2;
}

void initFolders(const std::string &output, Families &families) 
{
  std::string results = FileSystem::joinPaths(output, "results");
  FileSystem::mkdir(results, true);
  for (auto &family: families) {
    FileSystem::mkdir(FileSystem::joinPaths(results, family.name), true);
  }
}

void saveStats(const std::string &outputDir, double totalLibpllLL, double totalRecLL) 
{
  ParallelOfstream os(FileSystem::joinPaths(outputDir, "stats.txt"));
  os << "JointLL: " << totalLibpllLL + totalRecLL << std::endl;
  os << "LibpllLL: " << totalLibpllLL << std::endl;
  os << "RecLL: " << totalRecLL;
}

void optimizeStep(GeneRaxArguments &arguments, 
    RecModel recModel,
    Families &families,
    DTLRatesVector &rates,
    int sprRadius,
    int currentIteration,
    double &totalLibpllLL,
    double &totalRecLL,
    long &sumElapsedRates,
    long &sumElapsedSPR)
{
  Routines::optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, families, arguments.perSpeciesDTLRates, rates, sumElapsedRates);
  GeneTreeSearchMaster::optimizeGeneTrees(families, recModel, rates, arguments.output, "results",
      arguments.execPath, arguments.speciesTree, arguments.reconciliationOpt, arguments.rootedGeneTree,
      arguments.recWeight, true, true, sprRadius, currentIteration, useSplitImplem(), sumElapsedSPR);
  Routines::gatherLikelihoods(families, totalLibpllLL, totalRecLL);
}


void search(const Families &initialFamilies,
    GeneRaxArguments &arguments)

{
  Logger::timed << "Start SPR search" << std::endl;
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  long sumElapsedRates = 0;
  long sumElapsedSPR = 0;
  long sumElapsedLibpll = 0;
  DTLRatesVector rates(DTLRates(arguments.dupRate, arguments.lossRate, arguments.transferRate));
  Families currentFamilies = initialFamilies;
  int iteration = 0;
  bool randoms = Routines::createRandomTrees(arguments.output, currentFamilies); 
  if (randoms) {
    RaxmlMaster::runRaxmlOptimization(currentFamilies, arguments.output, arguments.execPath, iteration++, useSplitImplem(), sumElapsedLibpll);
    Routines::gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  }
  RecModel recModel = Arguments::strToRecModel(arguments.reconciliationModelStr);
  for (int i = 1; i <= arguments.maxSPRRadius; ++i) {
      optimizeStep(arguments, recModel, currentFamilies, rates, i, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  }
  saveStats(arguments.output, totalLibpllLL, totalRecLL);
  Routines::inferReconciliation(arguments.speciesTree, currentFamilies, recModel, rates, arguments.output);
  if (sumElapsedLibpll) {
    Logger::info << "Initial time spent on optimizing random trees: " << sumElapsedLibpll << "s" << std::endl;
  }
  Logger::info << "Time spent on optimizing rates: " << sumElapsedRates << "s" << std::endl;
  Logger::info << "Time spent on optimizing gene trees: " << sumElapsedSPR << "s" << std::endl;
  Logger::timed << "End of GeneRax execution" << std::endl;
}

void eval(const Families &initialFamilies,
    GeneRaxArguments &arguments)
{
  long dummy = 0;
  DTLRatesVector rates(DTLRates(arguments.dupRate, arguments.lossRate, arguments.transferRate));
  Families families = initialFamilies;
  RecModel recModel;
  bool randoms = Routines::createRandomTrees(arguments.output, families);
  if (randoms) {
    Logger::error << "You are running GeneRax in EVAL mode, but at least one starting gene tree was not provided. Aborting..." << std::endl;
    return;
  }
  recModel = Arguments::strToRecModel(arguments.reconciliationModelStr);
  Routines::optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, families, arguments.perSpeciesDTLRates, rates, dummy);
  int sprRadius = 0;
  int currentIteration = 0;
  GeneTreeSearchMaster::optimizeGeneTrees(families, recModel, rates, arguments.output, "results", arguments.execPath,
      arguments.speciesTree, arguments.reconciliationOpt, arguments.rootedGeneTree, arguments.recWeight,
      true, true, sprRadius, currentIteration, useSplitImplem(), dummy);
  
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  Routines::gatherLikelihoods(families, totalLibpllLL, totalRecLL);
  saveStats(arguments.output, totalLibpllLL, totalRecLL);
}


int generax_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm); 
  Logger::init();
  GeneRaxArguments arguments(argc, argv);
  ParallelContext::barrier();
  Logger::timed << "All cores started" << std::endl;
  srand(arguments.seed);
  FileSystem::mkdir(arguments.output, true);
  Logger::initFileOutput(FileSystem::joinPaths(arguments.output, "generax"));
  arguments.printCommand();
  arguments.printSummary();
  std::string labelledSpeciesTree = arguments.speciesTree + std::string(".labelled");
  LibpllParsers::labelRootedTree(arguments.speciesTree, labelledSpeciesTree);
  arguments.speciesTree = labelledSpeciesTree;
  Families initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::timed << "Number of gene families: " << initialFamilies.size() << std::endl;
  initFolders(arguments.output, initialFamilies);
  switch (arguments.strategy) {
  case SPR:
    search(initialFamilies, arguments);
    break;
  case EVAL:
    eval(initialFamilies, arguments);
    break;
  }
  Logger::close();
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
    return generax_main(argc, argv, comm);
  }
}

int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

