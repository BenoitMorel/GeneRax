#include "GeneRaxArguments.hpp"
#include <ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <limits>
#include <PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/DTLRates.hpp>
#include <treeSearch/JointTree.hpp>
#include <treeSearch/SPRSearch.hpp>
#include <IO/FileSystem.hpp>

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
    GeneRaxArguments &arguments) 
{
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
        arguments.check,
        false,
        rates.rates[0],
        rates.rates[1],
        rates.rates[2]);
    jointTree->optimizeParameters(true, false); // only optimize felsenstein likelihood
    double bestLoglk = jointTree->computeJointLoglk();
    jointTree->printLoglk();
    Logger::info << "Initial ll = " << bestLoglk << endl;
    while(SPRSearch::applySPRRound(*jointTree, 1, bestLoglk)) {} 
    Logger::info << "Final ll = " << bestLoglk << endl;

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
  

  // update DTL rates
  DTLRates rates(arguments.dupRate, arguments.lossRate, arguments.transferRate);
  if (!arguments.userDTLRates) {
    PerCoreGeneTrees geneTrees(initialFamilies);
    pll_rtree_t *speciesTree = LibpllEvaluation::readRootedFromFile(arguments.speciesTree); 
    rates = DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, arguments.reconciliationModel);
    pll_rtree_destroy(speciesTree, 0);
  }
  // update gene trees
  optimizeGeneTrees(initialFamilies, rates, arguments);


  Logger::timed << "End of GeneRax execution" << endl;
  ParallelContext::finalize();
  return 0;
}


#ifdef JOINTSEARCH_BUILD_AS_LIB

extern "C" int dll_main(int argc, char** argv, void* comm)
{
  return internal_main(argc, argv, comm);
}

#else

int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

#endif

