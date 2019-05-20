#include "GeneRaxArguments.hpp"
#include <ParallelContext.hpp>
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
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <sstream>



bool useSplitImplem() {
  return ParallelContext::getSize() > 4;
}

void schedule(const std::string &outputDir, const std::string &commandFile, bool splitImplem)
{
  std::vector<char *> argv;
  std::string exec = "mpi-scheduler";
  std::string implem = splitImplem ? "--split-scheduler" : "--fork-scheduler";
  
  std::string called_library = splitImplem ? "--static_scheduled_main" :  "/home/morelbt/github/GeneRax/build/bin/generaxslaves";
  std::string jobFailureFatal = "1";
  std::string threadsArg;
  std::string outputLogs = FileSystem::joinPaths(outputDir, "logs.txt");
  argv.push_back(const_cast<char *>(exec.c_str()));
  argv.push_back(const_cast<char *>(implem.c_str()));
  argv.push_back(const_cast<char *>(called_library.c_str()));
  argv.push_back(const_cast<char *>(commandFile.c_str()));
  argv.push_back(const_cast<char *>(outputDir.c_str()));
  argv.push_back(const_cast<char *>(jobFailureFatal.c_str()));
  argv.push_back(const_cast<char *>(threadsArg.c_str()));
  argv.push_back(const_cast<char *>(outputLogs.c_str()));
  MPI_Comm comm = MPI_COMM_WORLD;
  ParallelContext::barrier(); 
  if (splitImplem || ParallelContext::getRank() == 0) {
    mpi_scheduler_main(static_cast<int>(argv.size()), &argv[0], static_cast<void*>(&comm));
  }
  ParallelContext::barrier(); 
}

void raxmlMain(std::vector<FamiliesFileParser::FamilyInfo> &families,
    GeneRaxArguments &arguments,
    int iteration,
    long &sumElapsed)

{
  bool splitImplem = useSplitImplem();
  auto start = Logger::getElapsedSec();
  Logger::timed << "Starting raxml light step" << std::endl;
  std::stringstream outputDirName;
  outputDirName << "raxml_light_" << iteration;
  std::string outputDir = FileSystem::joinPaths(arguments.output, outputDirName.str());
  FileSystem::mkdir(outputDir, true);
  std::string commandFile = FileSystem::joinPaths(outputDir, "raxml_light_command.txt");
  auto geneTreeSizes = LibpllParsers::parallelGetTreeSizes(families);
  ParallelOfstream os(commandFile);
  for (size_t i = 0; i < families.size(); ++i) {
    auto &family = families[i];
    std::string familyOutput = FileSystem::joinPaths(arguments.output, "results");
    familyOutput = FileSystem::joinPaths(familyOutput, family.name);
    std::string geneTreePath = FileSystem::joinPaths(familyOutput, "geneTree.newick");
    std::string libpllModelPath = FileSystem::joinPaths(familyOutput, "libpllModel.txt");
    std::string outputStats = FileSystem::joinPaths(familyOutput, "raxml_light_stats.txt");
    auto taxa = geneTreeSizes[i];
    os << family.name << " ";
    os << 1 << " "; // cores
    os << taxa << " " ; // cost
    os << "raxmlLight" << " ";
    os << family.startingGeneTree << " ";
    os << family.alignmentFile << " ";
    os << family.libpllModel  << " ";
    os << geneTreePath << " ";
    os << libpllModelPath << " ";
    os << outputStats <<  std::endl;
    family.startingGeneTree = geneTreePath;
    family.statsFile = outputStats;
    family.libpllModel = libpllModelPath;
  }    
  os.close();
  schedule(outputDir, commandFile, splitImplem); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "End of raxml light step (after " << elapsed << "s)"  << std::endl;
}


void optimizeGeneTrees(std::vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates,
    GeneRaxArguments &arguments,
    RecModel recModel,
    bool enableRec,
    int sprRadius,
    int iteration,
    long &sumElapsed) 
{
  bool splitImplem = useSplitImplem();
  auto start = Logger::getElapsedSec();
  Logger::timed << "Starting SPR rounds with radius " << sprRadius << std::endl;
  std::stringstream outputDirName;
  outputDirName << "gene_optimization_" << iteration;
  std::string outputDir = FileSystem::joinPaths(arguments.output, outputDirName.str());
  FileSystem::mkdir(outputDir, true);
  std::string commandFile = FileSystem::joinPaths(outputDir, "opt_genes_command.txt");
  auto geneTreeSizes = LibpllParsers::parallelGetTreeSizes(families);
  ParallelOfstream os(commandFile);
  for (size_t i = 0; i < families.size(); ++i) {
    auto &family = families[i];
    std::string familyOutput = FileSystem::joinPaths(arguments.output, "results");
    familyOutput = FileSystem::joinPaths(familyOutput, family.name);
    std::string geneTreePath = FileSystem::joinPaths(familyOutput, "geneTree.newick");
    std::string outputStats = FileSystem::joinPaths(familyOutput, "stats.txt");
    auto taxa = geneTreeSizes[i];
    unsigned int cores = 1;
    if (sprRadius == 1) {
      cores = taxa / 2;;
    } else if (sprRadius == 2) {
      cores = taxa / 2;
    } else if (sprRadius >= 3) {
      cores = taxa;
    }
    if (cores == 0) {
      cores = 1;
    }
    if (!splitImplem) {
      cores = 1;
    }
    os << family.name << " ";
    os << cores << " "; // cores
    os << taxa << " " ; // cost
    os << "optimizeGeneTrees" << " ";
    os << family.startingGeneTree << " ";
    os << family.mappingFile << " ";
    os << family.alignmentFile << " ";
    os << arguments.speciesTree << " ";
    os << family.libpllModel  << " ";
    os << static_cast<int>(recModel)  << " ";
    os << static_cast<int>(arguments.reconciliationOpt)  << " ";
    os << static_cast<int>(arguments.rootedGeneTree)  << " ";
    os << arguments.recWeight  << " ";
    os << rates.rates[0]  << " ";
    os << rates.rates[1]  << " ";
    os << rates.rates[2]  << " ";
    os << static_cast<int>(enableRec)  << " ";
    os << sprRadius  << " ";
    os << geneTreePath << " ";
    os << outputStats <<  std::endl;
    family.startingGeneTree = geneTreePath;
    family.statsFile = outputStats;
  } 
  os.close();
  schedule(outputDir, commandFile, splitImplem); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "End of SPR rounds (after " << elapsed << "s)"  << std::endl;
}

void optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    RecModel recModel,
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates,
    long &sumElapsed) 
{
  if (userDTLRates) {
    return;
  }
  auto start = Logger::getElapsedSec();
  Logger::timed << "Start optimizing rates..." << std::endl;
  PerCoreGeneTrees geneTrees(families);
  pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
  rates = DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, recModel);
  //auto ratesVector = DTLOptimizer::optimizeDTLRatesVector(geneTrees, speciesTree, recModel);
  //Logger::info << "Optimized rates vector: " << ratesVector << std::endl;
  //Logger::info << "Optimized rates vector: " << ratesVector << std::endl;
  pll_rtree_destroy(speciesTree, 0);
  ParallelContext::barrier(); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "Finished optimizing rates: "
    << "D=" << rates.rates[0] << ", "
    << "L=" << rates.rates[1] << ", "
    << "T=" << rates.rates[2] << ", "
    << "Loglk=" << rates.ll 
    << " (after " << elapsed << "s)" << std::endl;
}

void inferReconciliation(
    const std::string &speciesTreeFile,
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    RecModel model,
    DTLRates &rates,
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
    evaluation.setRates(rates.rates[0], rates.rates[1], rates.rates[2]);
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
  /*
  for (unsigned int i = 0; i < speciesNodesCount; ++i) {
    Logger::info << dup_count[i] << std::endl;
    speciesTree->nodes[i]->length = dup_count[i] / static_cast<double>(families.size());
  }
  */
  auto speciesTreeStr = pll_rtree_export_newick(speciesTree->root, 0);
  Logger::info << "Species tree with average number of duplications (per familiy)" << std::endl;
  Logger::info << speciesTreeStr << std::endl;
  free(speciesTreeStr);
  pll_rtree_destroy(speciesTree, 0);
}


RecModel testRecModel(const std::string &speciesTreeFile,
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates)

{
  Logger::info << "Testing for transfers..." << std::endl;
  DTLRates transferRates;
  DTLRates noTransferRates;
  long testTime = 0;
  optimizeRates(false, speciesTreeFile, UndatedDL, families, noTransferRates, testTime);
  optimizeRates(false, speciesTreeFile, UndatedDTL, families, transferRates, testTime);
  double ll_dtl = transferRates.ll;
  double ll_dl = noTransferRates.ll;

  double aic_dtl = 2.0 - 2.0 * ll_dtl;
  double aic_dl = -2.0 * ll_dl;
  double bic_dtl = log(double(families.size())) - 2.0 * ll_dtl;
  double bic_dl = -2.0 * ll_dl;
  Logger::info << "AIC dl: " << aic_dl << std::endl;
  Logger::info << "AIC dtl: " << aic_dtl << std::endl;
  Logger::info << "BIC dl: " << bic_dl << std::endl;
  Logger::info << "BIC dtl: " << bic_dtl << std::endl;
  if (bic_dtl < bic_dl) {
    Logger::info << "According to BIC score, there were transfers between analyzed species." << std::endl;
    rates = transferRates;
    return UndatedDTL;
  } else {
    Logger::info << "According to BIC score, there were NO transfers between analyzed species." << std::endl;
    rates = noTransferRates;
    return UndatedDL;
  }
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
    DTLRates &rates,
    int sprRadius,
    int currentIteration,
    double &totalLibpllLL,
    double &totalRecLL,
    long &sumElapsedRates,
    long &sumElapsedSPR)
{
  optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, families, rates, sumElapsedRates);
  optimizeGeneTrees(families, rates, arguments, recModel, true, sprRadius, currentIteration, sumElapsedSPR);
  gatherLikelihoods(families, totalLibpllLL, totalRecLL);
}


void search(const std::vector<FamiliesFileParser::FamilyInfo> &initialFamilies,
    GeneRaxArguments &arguments)

{
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  long sumElapsedRates = 0;
  long sumElapsedSPR = 0;
  long sumElapsedLibpll = 0;

  DTLRates rates(arguments.dupRate, arguments.lossRate, arguments.transferRate);
  std::vector<FamiliesFileParser::FamilyInfo> currentFamilies = initialFamilies;

  bool randoms = createRandomTrees(arguments.output, currentFamilies); 
  int iteration = 0;
  if (randoms) {
    raxmlMain(currentFamilies, arguments, iteration++, sumElapsedLibpll);
    gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  }
  bool autoDetectRecModel;
  RecModel recModel;
  if (arguments.reconciliationModelStr == "AutoDetect") {
    autoDetectRecModel = true;
    recModel = UndatedDL;
  } else {
    autoDetectRecModel = false;
    recModel = Arguments::strToRecModel(arguments.reconciliationModelStr);
  }
  
  optimizeStep(arguments, recModel, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  optimizeStep(arguments, recModel, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  

  optimizeStep(arguments, recModel, currentFamilies, rates, 2, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  optimizeStep(arguments, recModel, currentFamilies, rates, 3, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);

  if (autoDetectRecModel) {
    recModel = testRecModel(arguments.speciesTree, currentFamilies, rates);
    optimizeStep(arguments, recModel, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
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
  DTLRates rates(arguments.dupRate, arguments.lossRate, arguments.transferRate);
  std::vector<FamiliesFileParser::FamilyInfo> families = initialFamilies;
  bool autoDetectRecModel;
  RecModel recModel;
  if (arguments.reconciliationModelStr == "AutoDetect") {
    autoDetectRecModel = true;
    recModel = UndatedDL;
  } else {
    autoDetectRecModel = false;
    recModel = Arguments::strToRecModel(arguments.reconciliationModelStr);
  }
  if (!autoDetectRecModel) {
    optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, families, rates, dummy);
  } else {
    recModel = testRecModel(arguments.speciesTree, families, rates);
  }
  int sprRadius = 0;
  int currentIteration = 0;
  optimizeGeneTrees(families, rates, arguments, recModel, true, sprRadius, currentIteration, dummy);
  
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  gatherLikelihoods(families, totalLibpllLL, totalRecLL);
  saveStats(arguments.output, totalLibpllLL, totalRecLL);
}


int internal_main(int argc, char** argv, void* comm)
{
  // the order is very important
  ParallelContext::init(comm); 
  Logger::init();
  GeneRaxArguments arguments(argc, argv);
  FileSystem::mkdir(arguments.output, true);
  Logger::initFileOutput(FileSystem::joinPaths(arguments.output, "generax"));
  
  arguments.printCommand();
  arguments.printSummary();

  std::vector<FamiliesFileParser::FamilyInfo> initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << std::endl;
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



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

