#include "GeneRaxArguments.hpp"
#include <ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <limits>
#include <PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/DTLRates.hpp>
#include <treeSearch/JointTree.hpp>
#include <treeSearch/SPRSearch.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <sstream>

using namespace std;


void optimizeGeneTrees(vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates,
    GeneRaxArguments &arguments,
    bool enableRec,
    int sprRadius,
    int iteration) 
{
  stringstream outputDirName;
  outputDirName << "gene_optimization_" << iteration;
  string commandFile = "/home/morelbt/github/phd_experiments/command.txt";
  string outputDir = FileSystem::joinPaths(arguments.output, outputDirName.str());
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
    } else if (sprRadius == 3) {
      cores = taxa;
    }
    cores = max(1, cores);
    //int cores = 8;// todobenoit use some smarter formula later
    //Logger::info << taxa << " " << cores << endl;
    os << family.name << " ";
    os << cores << " "; // cores
    os << taxa << " " ; // cost
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
  
  vector<char *> argv;
  string exec = "mpi-scheduler";
  string implem = "--split-scheduler" ;
  string library = "--static_scheduled_main"; //"/home/morelbt/github/GeneRax/build/src/generax/libgenerax_optimize_gene_trees.so";
  string jobFailureFatal = "1";
  string threadsArg;
  string outputLogs = FileSystem::joinPaths(outputDir, "logs.txt");
  FileSystem::mkdir(outputDir, true);
  argv.push_back((char *)exec.c_str());
  argv.push_back((char *)implem.c_str());
  argv.push_back((char *)library.c_str());
  argv.push_back((char *)commandFile.c_str());
  argv.push_back((char *)outputDir.c_str());
  argv.push_back((char *)jobFailureFatal.c_str());
  argv.push_back((char *)threadsArg.c_str());
  argv.push_back((char *)outputLogs.c_str());
  MPI_Comm comm = MPI_COMM_WORLD;
  Logger::timed << "Starting SPR rounds with radius " << sprRadius << endl;
  ParallelContext::barrier(); 
  mpi_scheduler_main(argv.size(), &argv[0], (void*)&comm);
}

void optimizeRates(bool userDTLRates, 
    const string &speciesTreeFile,
    RecModel recModel,
    vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates) 
{
  if (!userDTLRates) {
    PerCoreGeneTrees geneTrees(families);
    pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
    rates = DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, recModel);
    pll_rtree_destroy(speciesTree, 0);
    ParallelContext::barrier(); 
  }
}

RecModel testRecModel(const string &speciesTreeFile,
    vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates)

{
  DTLRates transferRates;
  DTLRates noTransferRates;
  
  optimizeRates(false, speciesTreeFile, UndatedDL, families, noTransferRates);
  optimizeRates(false, speciesTreeFile, UndatedDTL, families, transferRates);
  double myScoreTransfers = transferRates.ll / families.size();
  double myScoreNoTransfers = noTransferRates.ll / families.size();
  double ratio = fabs(myScoreNoTransfers / myScoreTransfers);
  Logger::info << "No transfer: " << noTransferRates.ll << " " << myScoreNoTransfers << endl;
  Logger::info << "With transfer: " << transferRates.ll << " " << myScoreTransfers << endl;
  Logger::info << "Ratio: " << ratio << endl;
  if (ratio > 1.1) {  
    Logger::info << "I think there are transfers" << endl;
    rates = transferRates;
    return UndatedDTL;
  } else {
    Logger::info << "I think there are NO transfers" << endl;
    rates = noTransferRates;
    return UndatedDL;
  }
}


void gatherLikelihoods(vector<FamiliesFileParser::FamilyInfo> &families,
    double &totalLibpllLL,
    double &totalRecLL)
{
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
  
  double totalLibpllLL = 0.0;
  double totalRecLL = 0.0;
  
  DTLRates rates(arguments.dupRate, arguments.lossRate, arguments.transferRate);
  vector<FamiliesFileParser::FamilyInfo> currentFamilies = initialFamilies;
  bool randoms = createRandomTrees(arguments.output, currentFamilies); 
  int iteration = 0;
  if (randoms) {
    optimizeGeneTrees(currentFamilies, rates, arguments, false, 1, iteration++);
  }
  RecModel recModel = arguments.reconciliationModel;

  if (arguments.autodetectDTLModel) {
    recModel = testRecModel(arguments.speciesTree, currentFamilies, rates);
    arguments.reconciliationModel = recModel;
  } else {
    optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, currentFamilies, rates);
  }

  optimizeGeneTrees(currentFamilies, rates, arguments, true, 1, iteration++);
  gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, true, 1, iteration++);
  gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, true, 2, iteration++);
  gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  //optimizeRates(arguments.userDTLRates, arguments.speciesTree, recModel, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, true, 3, iteration++);
  gatherLikelihoods(currentFamilies, totalLibpllLL, totalRecLL);
  saveStats(arguments.output, totalLibpllLL, totalRecLL);
  Logger::timed << "End of GeneRax execution" << endl;
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

