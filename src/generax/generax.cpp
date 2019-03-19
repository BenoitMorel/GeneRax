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
    string geneTreePath = FileSystem::joinPaths(arguments.output, family.name);
    geneTreePath = FileSystem::joinPaths(geneTreePath, "geneTree.newick");
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
    os << sprRadius  << " ";
    os << geneTreePath << endl;
    family.startingGeneTree = geneTreePath;
  } 
  os.close();
  
  vector<char *> argv;
  string exec = "mpi-scheduler";
  string implem = "--split-scheduler" ;
  string library = "--static_scheduled_main"; //"/home/morelbt/github/GeneRax/build/src/generax/libgenerax_optimize_gene_trees.so";
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
  mpi_scheduler_main(argv.size(), &argv[0], (void*)&comm);
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

void createRandomTrees(const string &geneRaxOutputDir, vector<FamiliesFileParser::FamilyInfo> &families)
{
  string startingTreesDir = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
  bool startingTreesDirCreated = false;
  for (auto &family: families) {
    if (family.startingGeneTree == "__random__") {
        if (!startingTreesDirCreated) {
          FileSystem::mkdir(startingTreesDir, true);
          startingTreesDir = true;
        } 
        family.startingGeneTree = FileSystem::joinPaths(geneRaxOutputDir, family.name + ".newick");
        if (ParallelContext::getRank() == 0) {
          LibpllEvaluation::createAndSaveRandomTree(family.alignmentFile, family.libpllModel, family.startingGeneTree);
        }
    }
  }
  ParallelContext::barrier();
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
  createRandomTrees(arguments.output, currentFamilies); 
  int iteration = 0; 
  optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 1, iteration++);
  optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 1, iteration++);
  optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 2, iteration++);
  //optimizeRates(arguments, currentFamilies, rates);
  optimizeGeneTrees(currentFamilies, rates, arguments, 3, iteration++);
  //optimizeRates(arguments, currentFamilies, rates);
  Logger::timed << "End of GeneRax execution" << endl;
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

