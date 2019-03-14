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
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>

using namespace std;

void getTreeStrings(const string &filename, vector<string> &treeStrings) 
{
  string geneTreeString;
  if (filename == "__random__" || filename.size() == 0) {
    treeStrings.push_back("__random__");
    return;
  }
  ifstream treeStream(filename);
  while(getline(treeStream, geneTreeString)) {
    geneTreeString.erase(remove(geneTreeString.begin(), geneTreeString.end(), '\n'), geneTreeString.end());
    if (geneTreeString.empty()) {
      continue;
    }
    treeStrings.push_back(geneTreeString);
  }
}

void optimizeGeneTrees(vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates,
    GeneRaxArguments &arguments,
    int sprRadius) 
{

/*
  vector<char *> argv;
  string exec = "mpi-scheduler";
  string implem = "--split-scheduler" ;
  string library = "/home/morelbt/github/GeneRax/build/src/generax/libgenerax_optimize_gene_trees.so";
  string commandFile = "/home/morelbt/github/phd_experiments/command.txt";
  string outputDir = FileSystem::joinPaths(arguments.output, "scheduler_run");
  string jobFailureFatal = "1";
  FileSystem::mkdir(outputDir, true);
  argv.push_back((char *)exec.c_str());
  argv.push_back((char *)implem.c_str());
  argv.push_back((char *)library.c_str());
  argv.push_back((char *)commandFile.c_str());
  argv.push_back((char *)outputDir.c_str());
  argv.push_back((char *)jobFailureFatal.c_str());
  MPI_Comm comm = MPI_COMM_WORLD;
  ParallelContext::barrier(); 
  cerr << "Begin rank " << ParallelContext::getRank() << endl;
  mpi_scheduler_main(argv.size(), &argv[0], (void*)&comm);
  */
  Logger::info << "Optimize gene trees with rates " << rates << endl;
  double totalInitialLL = 0.0;
  double totalFinalLL = 0.0;
  for (auto &family: families) {
    Logger::info << "Treating " << family.name << endl;
    vector<string> geneTreeStrings;
    getTreeStrings(family.startingGeneTree, geneTreeStrings);
    assert(geneTreeStrings.size() == 1);
    auto jointTree = make_shared<JointTree>(geneTreeStrings[0],
        family.alignmentFile,
        arguments.speciesTree,
        family.mappingFile,
        arguments.libpllModel,
        arguments.reconciliationModel,
        arguments.reconciliationOpt,
        arguments.rootedGeneTree,
        false,
        false,
        rates.rates[0],
        rates.rates[1],
        rates.rates[2]);
    jointTree->optimizeParameters(true, false); // only optimize felsenstein likelihood
    double bestLoglk = jointTree->computeJointLoglk();
    totalInitialLL += bestLoglk;
    jointTree->printLoglk();
    Logger::info << "Initial ll = " << bestLoglk << endl;
    while(SPRSearch::applySPRRound(*jointTree, sprRadius, bestLoglk)) {} 
    totalFinalLL += bestLoglk;
    Logger::info << "Final ll = " << bestLoglk << endl;
    string geneTreePath = FileSystem::joinPaths(arguments.output, family.name);
    geneTreePath = FileSystem::joinPaths(geneTreePath, "geneTree.newick");
    jointTree->save(geneTreePath, false);
    family.startingGeneTree = geneTreePath;
  }
  Logger::info << "Total initial and final ll: " << totalInitialLL << " " << totalFinalLL << endl;
  cerr << "End rank " << ParallelContext::getRank() << endl;
  ParallelContext::barrier(); 
}

void optimizeRates(const GeneRaxArguments &arguments,
    vector<FamiliesFileParser::FamilyInfo> &families,
    DTLRates &rates) 
{
  if (!arguments.userDTLRates) {
    PerCoreGeneTrees geneTrees(families);
    pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(arguments.speciesTree); 
    rates = DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, arguments.reconciliationModel);
    pll_rtree_destroy(speciesTree, 0);
    ParallelContext::barrier(); 
  }
}

void initFolders(const string &output, vector<FamiliesFileParser::FamilyInfo> &families) 
{
  FileSystem::mkdir(output, true);
  for (auto &family: families) {
    FileSystem::mkdir(FileSystem::joinPaths(output, family.name), true);
  }
}

int internal_main(int argc, char** argv, void* comm)
{
  // the order is very important
  ParallelContext::init(comm); 
  Logger::init();
  GeneRaxArguments arguments(argc, argv);
  Logger::initFileOutput(arguments.output);
  
  arguments.printCommand();
  arguments.printSummary();

  vector<FamiliesFileParser::FamilyInfo> initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << endl;
  initFolders(arguments.output, initialFamilies);
  
  
  DTLRates rates(arguments.dupRate, arguments.lossRate, arguments.transferRate);
  vector<FamiliesFileParser::FamilyInfo> currentFamilies = initialFamilies;
  
  optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 1);
  optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 1);
  optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 2);
  optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 3);
  optimizeRates(arguments, currentFamilies, rates);
  Logger::timed << "End of GeneRax execution" << endl;
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}


