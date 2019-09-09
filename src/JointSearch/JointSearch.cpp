#include "JSArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <trees/JointTree.hpp>
#include <search/SPRSearch.hpp>
#include <IO/Logger.hpp>
#include <util/Scenario.hpp>
#include <maths/DTLRates.hpp>
#include <algorithm>
#include <limits>
#include <cmath>



/**
 *  @param: filename: a file with one newick std::string per line
 *  @param treeStrings: std::vector of the newick strings, formated such that 
 *    JointSearch can read them
 */
void getTreeStrings(const std::string &filename, std::vector<std::string> &treeStrings) 
{
  std::string geneTreeString;
  if (filename == "__random__" || filename.size() == 0) {
    treeStrings.push_back("__random__");
    return;
  }
  std::ifstream treeStream(filename);
  while(getline(treeStream, geneTreeString)) {
    geneTreeString.erase(remove(geneTreeString.begin(), geneTreeString.end(), '\n'), geneTreeString.end());
    if (geneTreeString.empty()) {
      continue;
    }
    treeStrings.push_back(geneTreeString);
  }
}

/**
 *  JointSearch main function
 *  @param argc:
 *  @param argv:
 *  @param comm: the communicator that JointSearch is allowed to use
 *    for its parallel context (relevant when called as a lib by another program)
 */
int internal_main(int argc, char** argv, void* comm)
{
  // the order is very important
  ParallelContext::init(comm); 
  Logger::init();
  JSArguments arguments(argc, argv);
  Logger::initFileOutput(arguments.output);
  
  arguments.printCommand();
  arguments.printSummary();
  
  std::vector<std::string> geneTreeStrings;
  getTreeStrings(arguments.geneTree, geneTreeStrings);
  
  bool firstRun  = true; 
  double bestLL = std::numeric_limits<double>::lowest();
  
  std::string bestTreeFile = arguments.output + ".newick";
  std::string allTreesFile = arguments.output + "_all" + ".newick";
  std::string eventCountsFile = arguments.output + ".events";
  std::string treeWithEventsFile = arguments.output + "_withevents.nhx";
  std::string statsFile = arguments.output + ".stats";
  for (auto &geneTreeString: geneTreeStrings) {
    double dupRate = 1;
    double lossRate = 1;
    double transferRate = 1;
    if (arguments.userDTLRates) {
      dupRate = arguments.dupRate;
      lossRate = arguments.lossRate;
      transferRate = arguments.transferRate;
    }
    auto jointTree = std::make_shared<JointTree>(geneTreeString,
        arguments.alignment,
        arguments.speciesTree,
        arguments.geneSpeciesMap,
        arguments.libpllModel,
        arguments.reconciliationModel,
        arguments.reconciliationOpt,
        arguments.rootedGeneTree,
        false, // prune species tree
        1.0,
        arguments.check,
        true, // optimize DTL rates?
        DTLRatesVector(DTLRates(dupRate, lossRate, transferRate))
        );
    jointTree->printInfo();;
    jointTree->optimizeParameters();
    double initialRecLL = jointTree->computeReconciliationLoglk();
    double initialLibpllLL = jointTree->computeLibpllLoglk();
    Logger::timed << "Starting search..." << std::endl;
    if (arguments.strategy == SPR) {
      if (!geneTreeString.size() or geneTreeString == "__random__") {
        jointTree->enableReconciliation(false);
        SPRSearch::applySPRSearch(*jointTree);
        jointTree->enableReconciliation(true);
      }
      SPRSearch::applySPRSearch(*jointTree);
    } else if (arguments.strategy == EVAL) {
    }
    Logger::timed << "End of search" << std::endl;
    jointTree->printLoglk();
    Logger::info << "Final tree hash: " << jointTree->getUnrootedTreeHash() << std::endl;
    if (!ParallelContext::getRank()) {
      double ll = jointTree->computeJointLoglk();
      assert(!std::isnan(ll));
      if (ll >= bestLL) {
        bestLL = ll;
        jointTree->save(bestTreeFile, false);
        ParallelOfstream stats(statsFile);
        stats << "initial_ll " << initialRecLL + initialLibpllLL << std::endl;
        stats << "initial_llrec " << initialRecLL << std::endl;
        stats << "initial_lllibpll " << initialLibpllLL << std::endl;
        stats << "ll " << bestLL << std::endl;
        stats << "llrec " << jointTree->computeReconciliationLoglk() << std::endl;
        stats << "lllibpll " << jointTree->computeLibpllLoglk() << std::endl;
        stats << "D " << jointTree->getRatesVector().getRates(0).rates[0] << std::endl;
        stats << "L " << jointTree->getRatesVector().getRates(0).rates[1] << std::endl;
        stats << "T " << jointTree->getRatesVector().getRates(0).rates[2] << std::endl;
        stats << "hash " << jointTree->getUnrootedTreeHash() << std::endl;
        stats << " " << std::endl;
      }
      jointTree->save(allTreesFile, !firstRun);
      Scenario scenario;
      jointTree->inferMLScenario(scenario);
      Logger::info << std::endl;
      scenario.saveEventsCounts(eventCountsFile);
      scenario.saveReconciliation(treeWithEventsFile, NHX);
    }
    firstRun = false;
  }  
  Logger::info << "Best tree: " + bestTreeFile << std::endl;
  Logger::info << "Best tree with events: " + treeWithEventsFile << std::endl;
  Logger::timed << "End of JointSearch execution" << std::endl;
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

