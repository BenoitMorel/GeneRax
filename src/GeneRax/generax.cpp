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
#include <treeSearch/SPRSearch.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <sstream>

using namespace std;

void schedule(const string &outputDir, const string &commandFile)
{
  vector<char *> argv;
  string exec = "mpi-scheduler";
  string implem = "--split-scheduler" ;
  string library = "--static_scheduled_main"; //"/home/morelbt/github/GeneRax/build/src/generax/libgenerax_optimize_gene_trees.so";
  string jobFailureFatal = "1";
  string threadsArg;
  string outputLogs = FileSystem::joinPaths(outputDir, "logs.txt");
  argv.push_back((char *)exec.c_str());
  argv.push_back((char *)implem.c_str());
  argv.push_back((char *)library.c_str());
  argv.push_back((char *)commandFile.c_str());
  argv.push_back((char *)outputDir.c_str());
  argv.push_back((char *)jobFailureFatal.c_str());
  argv.push_back((char *)threadsArg.c_str());
  argv.push_back((char *)outputLogs.c_str());
  MPI_Comm comm = MPI_COMM_WORLD;
  ParallelContext::barrier(); 
  mpi_scheduler_main(argv.size(), &argv[0], (void*)&comm);
}

void raxmlMain(vector<FamiliesFileParser::FamilyInfo> &families,
    GeneRaxArguments &arguments,
    int iteration,
    long &sumElapsed)

{
  auto start = Logger::getElapsedSec();
  Logger::timed << "Starting raxml light step" << endl;
  
  stringstream outputDirName;
  outputDirName << "raxml_light_" << iteration;
  string outputDir = FileSystem::joinPaths(arguments.output, outputDirName.str());
  FileSystem::mkdir(outputDir, true);
  string commandFile = FileSystem::joinPaths(outputDir, "raxml_light_command.txt");
  vector<int> geneTreeSizes = LibpllParsers::parallelGetTreeSizes(families);
  ParallelOfstream os(commandFile);
  for (size_t i = 0; i < families.size(); ++i) {
    auto &family = families[i];
    string familyOutput = FileSystem::joinPaths(arguments.output, "results");
    familyOutput = FileSystem::joinPaths(familyOutput, family.name);
    string geneTreePath = FileSystem::joinPaths(familyOutput, "geneTree.newick");
    string libpllModelPath = FileSystem::joinPaths(familyOutput, "libpllModel.txt");
    string outputStats = FileSystem::joinPaths(familyOutput, "raxml_light_stats.txt");
    int taxa = geneTreeSizes[i];
    os << family.name << " ";
    os << 1 << " "; // cores
    os << taxa << " " ; // cost
    os << "raxmlLight" << " ";
    os << family.startingGeneTree << " ";
    os << family.alignmentFile << " ";
    os << family.libpllModel  << " ";
    os << geneTreePath << " ";
    os << libpllModelPath << " ";
    os << outputStats <<  endl;
    family.startingGeneTree = geneTreePath;
    family.statsFile = outputStats;
    family.libpllModel = libpllModelPath;
  }    
  os.close();
  schedule(outputDir, commandFile); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "End of raxml light step (after " << elapsed << "s)"  << endl;
}


void optimizeGeneTrees(vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates,
    GeneRaxArguments &arguments,
    bool enableRec,
    int sprRadius,
    int iteration,
    long &sumElapsed) 
{
  auto start = Logger::getElapsedSec();
  Logger::timed << "Starting SPR rounds with radius " << sprRadius << endl;
  stringstream outputDirName;
  outputDirName << "gene_optimization_" << iteration;
  string outputDir = FileSystem::joinPaths(arguments.output, outputDirName.str());
  FileSystem::mkdir(outputDir, true);
  string commandFile = FileSystem::joinPaths(outputDir, "opt_genes_command.txt");
  vector<int> geneTreeSizes = LibpllParsers::parallelGetTreeSizes(families);
  ParallelOfstream os(commandFile);
  for (size_t i = 0; i < families.size(); ++i) {
    auto &family = families[i];
    string familyOutput = FileSystem::joinPaths(arguments.output, "results");
    familyOutput = FileSystem::joinPaths(familyOutput, family.name);
    string geneTreePath = FileSystem::joinPaths(familyOutput, "geneTree.newick");
    string outputStats = FileSystem::joinPaths(familyOutput, "stats.txt");
    int taxa = geneTreeSizes[i];
    int cores = 1;
    if (sprRadius == 1) {
      cores = taxa / 2;;
    } else if (sprRadius == 2) {
      cores = taxa / 2;
    } else if (sprRadius >= 3) {
      cores = taxa;
    }
    cores = max(1, cores);
    os << family.name << " ";
    os << cores << " "; // cores
    os << taxa << " " ; // cost
    os << "optimizeGeneTrees" << " ";
    os << family.startingGeneTree << " ";
    os << family.mappingFile << " ";
    os << family.alignmentFile << " ";
    os << arguments.speciesTree << " ";
    os << family.libpllModel  << " ";
    os << (int)arguments.reconciliationModel  << " ";
    os << (int)arguments.reconciliationOpt  << " ";
    os << (int)arguments.rootedGeneTree  << " ";
    os << rates.rates[0]  << " ";
    os << rates.rates[1]  << " ";
    os << rates.rates[2]  << " ";
    os << int(enableRec)  << " ";
    os << sprRadius  << " ";
    os << geneTreePath << " ";
    os << outputStats <<  endl;
    family.startingGeneTree = geneTreePath;
    family.statsFile = outputStats;
  } 
  os.close();
  schedule(outputDir, commandFile); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "End of SPR rounds (after " << elapsed << "s)"  << endl;
}

void optimizeRates(bool userDTLRates, 
    const string &speciesTreeFile,
    RecModel recModel,
    vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates,
    long &sumElapsed) 
{
  auto start = Logger::getElapsedSec();
  Logger::timed << "Start optimizing rates..." << endl;
  if (!userDTLRates) {
    PerCoreGeneTrees geneTrees(families);
    pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
    rates = DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, recModel);
    pll_rtree_destroy(speciesTree, 0);
    ParallelContext::barrier(); 
  }
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "Finished optimizing rates: "
    << "D=" << rates.rates[0] << ", "
    << "L=" << rates.rates[1] << ", "
    << "T=" << rates.rates[2] << ", "
    << "Loglk=" << rates.ll 
    << " (after " << elapsed << "s)" << endl;
}

RecModel testRecModel(const string &speciesTreeFile,
    vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates)

{
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
  Logger::info << "AIC dl: " << aic_dl << endl;
  Logger::info << "AIC dtl: " << aic_dtl << endl;
  Logger::info << "BIC dl: " << bic_dl << endl;
  Logger::info << "BIC dtl: " << bic_dtl << endl;
  if (aic_dtl < aic_dl) {
    Logger::info << "AIC thinks there are some transfers" << endl;
  } else {
    Logger::info << "AIC thinks there are NO transfers" << endl;
  }
  if (bic_dtl < bic_dl) {
    Logger::info << "BIC thinks there are some transfers" << endl;
  } else {
    Logger::info << "BIC thinks there are NO transfers" << endl;
  }
  
  // todobenoit
  return UndatedDTL;
}


void gatherLikelihoods(vector<FamiliesFileParser::FamilyInfo> &families,
    double &totalLibpllLL,
    double &totalRecLL)
{
  Logger::info << "Start gathering likelihoods... " << endl;
  totalRecLL = 0.0;
  totalLibpllLL = 0.0;
  int familiesNumber = families.size();
  for (int i = ParallelContext::getBegin(familiesNumber); i < ParallelContext::getEnd(familiesNumber); ++i) {
    auto &family = families[i];
    ifstream is(family.statsFile);
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
  Logger::info << "rec = " << totalRecLL << endl;
}

void initFolders(const string &output, vector<FamiliesFileParser::FamilyInfo> &families) 
{
  string results = FileSystem::joinPaths(output, "results");
  FileSystem::mkdir(results, true);
  for (auto &family: families) {
    FileSystem::mkdir(FileSystem::joinPaths(results, family.name), true);
  }
}

bool createRandomTrees(const string &geneRaxOutputDir, vector<FamiliesFileParser::FamilyInfo> &families)
{
  string startingTreesDir = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
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

void saveStats(const string &outputDir, double totalLibpllLL, double totalRecLL) 
{
  ParallelOfstream os(FileSystem::joinPaths(outputDir, "stats.txt"));
  os << "JointLL: " << totalLibpllLL + totalRecLL << endl;
  os << "LibpllLL: " << totalLibpllLL << endl;
  os << "RecLL:" << totalRecLL;
}

void optimizeStep(GeneRaxArguments &arguments, 
    vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates,
    int sprRadius,
    int currentIteration,
    double &totalLibpllLL,
    double &totalRecLL,
    long &sumElapsedRates,
    long &sumElapsedSPR)
{
  optimizeRates(arguments.userDTLRates, arguments.speciesTree, arguments.reconciliationModel, families, rates, sumElapsedRates);
  optimizeGeneTrees(families, rates, arguments, true, sprRadius, currentIteration, sumElapsedSPR);
  gatherLikelihoods(families, totalLibpllLL, totalRecLL);
}


void search(const vector<FamiliesFileParser::FamilyInfo> &initialFamilies,
    GeneRaxArguments &arguments)

{
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  long sumElapsedRates = 0;
  long sumElapsedSPR = 0;
  long sumElapsedLibpll = 0;

  DTLRates rates(arguments.dupRate, arguments.lossRate, arguments.transferRate);
  vector<FamiliesFileParser::FamilyInfo> currentFamilies = initialFamilies;

  bool randoms = createRandomTrees(arguments.output, currentFamilies); 
  int iteration = 0;
  if (randoms) {
    raxmlMain(currentFamilies, arguments, iteration++, sumElapsedLibpll);
    gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  }
  RecModel initialRecModel = arguments.reconciliationModel;
  if (initialRecModel == UndatedDTL) {
    arguments.reconciliationModel = UndatedDL;
  }
 
  optimizeStep(arguments, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  optimizeStep(arguments, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  

  optimizeStep(arguments, currentFamilies, rates, 2, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  optimizeStep(arguments, currentFamilies, rates, 3, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);

  if (initialRecModel == UndatedDTL) {
    DTLRates rates;
    testRecModel(arguments.speciesTree, currentFamilies, rates);
    arguments.reconciliationModel = UndatedDTL;
    optimizeStep(arguments, currentFamilies, rates, 1, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);
  }
  
  optimizeStep(arguments, currentFamilies, rates, arguments.maxSPRRadius, iteration++, totalLibpllLL, totalRecLL, sumElapsedRates, sumElapsedSPR);

  saveStats(arguments.output, totalLibpllLL, totalRecLL);
  if (sumElapsedLibpll) {
    Logger::info << "Initial time spent on optimizing random trees: " << sumElapsedLibpll << "s" << endl;
  }
  Logger::info << "Time spent on optimizing rates: " << sumElapsedRates << "s" << endl;
  Logger::info << "Time spent on optimizing gene trees: " << sumElapsedSPR << "s" << endl;
  Logger::timed << "End of GeneRax execution" << endl;

}

void eval(const vector<FamiliesFileParser::FamilyInfo> &initialFamilies,
    GeneRaxArguments &arguments)
{
  long dummy = 0;
  DTLRates rates(arguments.dupRate, arguments.lossRate, arguments.transferRate);
  vector<FamiliesFileParser::FamilyInfo> families = initialFamilies;
  optimizeRates(arguments.userDTLRates, arguments.speciesTree, arguments.reconciliationModel, families, rates, dummy);
  int sprRadius = 0;
  int currentIteration = 0;
  optimizeGeneTrees(families, rates, arguments, true, sprRadius, currentIteration, dummy);
  
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  gatherLikelihoods(families, totalLibpllLL, totalRecLL);

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

  vector<FamiliesFileParser::FamilyInfo> initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << endl;
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

