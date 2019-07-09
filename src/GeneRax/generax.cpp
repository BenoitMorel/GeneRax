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
#include <trees/JointTree.hpp>
#include <search/SPRSearch.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <sstream>
#include "GeneRaxSlaves.hpp"
#include <parallelization/Scheduler.hpp>
#include <routines/RaxmlMaster.hpp>
#include <routines/GeneTreeSearchMaster.hpp>

bool useSplitImplem() {
  return ParallelContext::getSize() > 2;
}


void optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    RecModel recModel,
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    bool perSpeciesRates, 
    DTLRatesVector &rates,
    long &sumElapsed) 
{
  if (userDTLRates) {
    return;
  }
  auto start = Logger::getElapsedSec();
  Logger::timed << "Start optimizing rates..." << std::endl;
  PerCoreGeneTrees geneTrees(families);
  bool ok = geneTrees.checkMappings(speciesTreeFile); 
  if (!ok) {
    Logger::info << "INVALID MAPPINGS" << std::endl;
    ParallelContext::abort(1);
  }
  pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
  if (!perSpeciesRates) {
    rates = DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, recModel);
  } else {
    rates = DTLOptimizer::optimizeDTLRatesVector(geneTrees, speciesTree, recModel, &rates);
  }
  pll_rtree_destroy(speciesTree, 0);
  ParallelContext::barrier(); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "Finished optimizing rates: "
    << "Loglk=" << rates.getLL() 
    << " (after " << elapsed << "s)" << std::endl;
}

void inferReconciliation(
    const std::string &speciesTreeFile,
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    RecModel model,
    DTLRatesVector &rates,
    const std::string &outputDir
    )
{
  pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
  PerCoreGeneTrees geneTrees(families);
  std::string reconciliationsDir = FileSystem::joinPaths(outputDir, "reconciliations");
  FileSystem::mkdir(reconciliationsDir, true);
  auto speciesNodesCount = speciesTree->tip_count + speciesTree->inner_count;
  std::vector<double> dup_count(speciesNodesCount, 0.0);
  ParallelContext::barrier();
  for (auto &tree: geneTrees.getTrees()) {
    std::string eventCountsFile = FileSystem::joinPaths(reconciliationsDir, tree.name + "_eventCounts.txt");
    std::string treeWithEventsFile = FileSystem::joinPaths(reconciliationsDir, tree.name + "_reconciliated.nhx");
    Scenario scenario;
    ReconciliationEvaluation evaluation(speciesTree, tree.mapping, model, true);
    evaluation.setRates(rates);
    evaluation.evaluate(tree.tree);
    evaluation.inferMLScenario(scenario);
    scenario.saveEventsCounts(eventCountsFile, false);
    scenario.saveTreeWithEvents(treeWithEventsFile, false);
    for (auto &event: scenario.getEvents()) {
      if (event.type == Scenario::D) {
        dup_count[event.speciesNode]++;
      } 
    }
  }
  ParallelContext::sumVectorDouble(dup_count);
  auto initialSpeciesTreeStr = pll_rtree_export_newick(speciesTree->root, 0);
  free(initialSpeciesTreeStr);
  auto speciesTreeStr = pll_rtree_export_newick(speciesTree->root, 0);
  free(speciesTreeStr);
  pll_rtree_destroy(speciesTree, 0);
}



void gatherLikelihoods(std::vector<FamiliesFileParser::FamilyInfo> &families,
    double &totalLibpllLL,
    double &totalRecLL)
{
  Logger::info << "Start gathering likelihoods... " << std::endl;
  totalRecLL = 0.0;
  totalLibpllLL = 0.0;
  unsigned int familiesNumber = static_cast<unsigned int>(families.size());
  for (auto i = ParallelContext::getBegin(familiesNumber); i < ParallelContext::getEnd(familiesNumber); ++i) {
    auto &family = families[i];
    std::ifstream is(family.statsFile);
    double libpllLL = 0.0;
    double recLL = 0.0;
    is >> libpllLL;
    is >> recLL;
    totalRecLL += recLL;
    totalLibpllLL += libpllLL;
  }
  ParallelContext::sumDouble(totalRecLL);
  ParallelContext::sumDouble(totalLibpllLL);
  Logger::info << "Likelihoods: ";
  Logger::info << "joint = " << totalLibpllLL + totalRecLL << ", ";
  Logger::info << "libpll = " << totalLibpllLL << ", ";
  Logger::info << "rec = " << totalRecLL << std::endl;
}

void initFolders(const std::string &output, std::vector<FamiliesFileParser::FamilyInfo> &families) 
{
  std::string results = FileSystem::joinPaths(output, "results");
  FileSystem::mkdir(results, true);
  for (auto &family: families) {
    FileSystem::mkdir(FileSystem::joinPaths(results, family.name), true);
  }
}

bool createRandomTrees(const std::string &geneRaxOutputDir, std::vector<FamiliesFileParser::FamilyInfo> &families)
{
  std::string startingTreesDir = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
  bool startingTreesDirCreated = false;
  for (auto &family: families) {
    if (family.startingGeneTree == "__random__") {
        if (!startingTreesDirCreated) {
          FileSystem::mkdir(startingTreesDir, true);
          startingTreesDirCreated = true;
        } 
        family.startingGeneTree = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
        family.startingGeneTree = FileSystem::joinPaths(family.startingGeneTree, family.name + ".newick");
        if (ParallelContext::getRank() == 0) {
          LibpllEvaluation::createAndSaveRandomTree(family.alignmentFile, family.libpllModel, family.startingGeneTree);
        }
    }
  }
  ParallelContext::barrier();
  return startingTreesDirCreated;
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
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRatesVector &rates,
    int sprRadius,
    int currentIteration,
    double &totalLibpllLL,
    double &totalRecLL,
    long &sumElapsedRates,
    long &sumElapsedSPR)
{
  optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, families, arguments.perSpeciesDTLRates, rates, sumElapsedRates);
  GeneTreeSearchMaster::optimizeGeneTrees(families, 
      recModel,
      rates, 
      arguments.output,
      arguments.execPath,
      arguments.speciesTree,
      arguments.reconciliationOpt,
      arguments.rootedGeneTree,
      arguments.recWeight,
      true,
      sprRadius,
      currentIteration,
      useSplitImplem(),
      sumElapsedSPR);
  gatherLikelihoods(families, totalLibpllLL, totalRecLL);
}


void search(const std::vector<FamiliesFileParser::FamilyInfo> &initialFamilies,
    GeneRaxArguments &arguments)

{
  Logger::timed << "Start SPR search" << std::endl;
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  long sumElapsedRates = 0;
  long sumElapsedSPR = 0;
  long sumElapsedLibpll = 0;

  DTLRatesVector rates(DTLRates(arguments.dupRate, arguments.lossRate, arguments.transferRate));
  std::vector<FamiliesFileParser::FamilyInfo> currentFamilies = initialFamilies;

  bool randoms = createRandomTrees(arguments.output, currentFamilies); 
  int iteration = 0;
  
  if (randoms) {
    RaxmlMaster::runRaxmlOptimization(currentFamilies, arguments.output, arguments.execPath, iteration++, useSplitImplem(), sumElapsedLibpll);
    gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  }
  RecModel recModel = Arguments::strToRecModel(arguments.reconciliationModelStr);
  
  optimizeStep(arguments, recModel, currentFamilies, rates, 0, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  optimizeStep(arguments, recModel, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  optimizeStep(arguments, recModel, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  
  if (arguments.maxSPRRadius >= 2) {
    optimizeStep(arguments, recModel, currentFamilies, rates, 2, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  }
  if (arguments.maxSPRRadius >= 3) {
    optimizeStep(arguments, recModel, currentFamilies, rates, 3, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  }
  optimizeStep(arguments, recModel, currentFamilies, rates, arguments.maxSPRRadius, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);

  saveStats(arguments.output, totalLibpllLL, totalRecLL);
  inferReconciliation(arguments.speciesTree, currentFamilies, recModel, rates, arguments.output);
  if (sumElapsedLibpll) {
    Logger::info << "Initial time spent on optimizing random trees: " << sumElapsedLibpll << "s" << std::endl;
  }
  Logger::info << "Time spent on optimizing rates: " << sumElapsedRates << "s" << std::endl;
  Logger::info << "Time spent on optimizing gene trees: " << sumElapsedSPR << "s" << std::endl;
  Logger::timed << "End of GeneRax execution" << std::endl;

}

void eval(const std::vector<FamiliesFileParser::FamilyInfo> &initialFamilies,
    GeneRaxArguments &arguments)
{
  long dummy = 0;
  DTLRatesVector rates(DTLRates(arguments.dupRate, arguments.lossRate, arguments.transferRate));
  std::vector<FamiliesFileParser::FamilyInfo> families = initialFamilies;
  RecModel recModel;
  bool randoms = createRandomTrees(arguments.output, families);
  if (randoms) {
    Logger::info << "[Warning] You are running GeneRax in EVAL mode, but at least one starting gene tree was not provided (a random tree was generated instead)" << std::endl;
  }
  recModel = Arguments::strToRecModel(arguments.reconciliationModelStr);
  optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, families, arguments.perSpeciesDTLRates, rates, dummy);
  int sprRadius = 0;
  int currentIteration = 0;
  GeneTreeSearchMaster::optimizeGeneTrees(families, 
      recModel,
      rates, 
      arguments.output,
      arguments.execPath,
      arguments.speciesTree,
      arguments.reconciliationOpt,
      arguments.rootedGeneTree,
      arguments.recWeight,
      true,
      sprRadius,
      currentIteration,
      useSplitImplem(),
      dummy);
  
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  gatherLikelihoods(families, totalLibpllLL, totalRecLL);
  saveStats(arguments.output, totalLibpllLL, totalRecLL);
}


int generax_main(int argc, char** argv, void* comm)
{
  // the order is very important
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

  std::vector<FamiliesFileParser::FamilyInfo> initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
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

