#include <IO/Arguments.hpp>
#include <ParallelContext.hpp>
#include <treeSearch/JointTree.hpp>
#include <treeSearch/SPRSearch.hpp>
#include <IO/Logger.hpp>
#include <Scenario.hpp>

#include <algorithm>
#include <limits>

using namespace std;


/**
 *  @param: filename: a file with one newick string per line
 *  @param treeStrings: vector of the newick strings, formated such that 
 *    JointSearch can read them
 */
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
  Arguments arguments(argc, argv);
  Logger::initFileOutput(arguments.output);
  
  arguments.printCommand();
  arguments.printSummary();
  
  vector<string> geneTreeStrings;
  getTreeStrings(arguments.geneTree, geneTreeStrings);
  
  bool firstRun  = true; 
  double bestLL = numeric_limits<double>::lowest();
  
  string bestTreeFile = arguments.output + ".newick";
  string allTreesFile = arguments.output + "_all" + ".newick";
  string eventCountsFile = arguments.output + ".events";
  string treeWithEventsFile = arguments.output + "_withevents.nhx";
  string statsFile = arguments.output + ".stats";
  for (auto &geneTreeString: geneTreeStrings) {
    double dupRate = 1;
    double lossRate = 1;
    double transferRate = 1;
    if (arguments.userDTLRates) {
      dupRate = arguments.dupRate;
      lossRate = arguments.lossRate;
      transferRate = arguments.transferRate;
    }
    auto jointTree = make_shared<JointTree>(geneTreeString,
        arguments.alignment,
        arguments.speciesTree,
        arguments.geneSpeciesMap,
        arguments.libpllModel,
        arguments.reconciliationModel,
        arguments.reconciliationOpt,
        arguments.rootedGeneTree,
        arguments.check,
        dupRate,
        lossRate,
        transferRate);
    jointTree->printInfo();;
    jointTree->optimizeParameters();
    double initialRecLL = jointTree->computeReconciliationLoglk();
    double initialLibpllLL = jointTree->computeLibpllLoglk();
    Logger::timed << "Starting search..." << endl;
    if (arguments.strategy == "SPR") {
      if (!geneTreeString.size() or geneTreeString == "__random__") {
        Logger::info << "plip" << endl;
        jointTree->enableReconciliation(false);
        SPRSearch::applySPRSearch(*jointTree);
        Logger::info << "plop" << endl;
        jointTree->enableReconciliation(true);
      }
      SPRSearch::applySPRSearch(*jointTree);
    } else if (arguments.strategy == "EVAL") {
    }
    Logger::timed << "End of search" << endl;
    jointTree->printLoglk();
    Logger::info << "Final tree hash: " << jointTree->getUnrootedTreeHash() << endl;
    if (!ParallelContext::getRank()) {
      double ll = jointTree->computeJointLoglk();
      assert(!isnan(ll));
      if (ll >= bestLL) {
        bestLL = ll;
        jointTree->save(bestTreeFile, false);
        ParallelOfstream stats(statsFile);
        stats << "initial_ll " << initialRecLL + initialLibpllLL << endl;
        stats << "initial_llrec " << initialRecLL << endl;
        stats << "initial_lllibpll " << initialLibpllLL << endl;
        stats << "ll " << bestLL << endl;
        stats << "llrec " << jointTree->computeReconciliationLoglk() << endl;
        stats << "lllibpll " << jointTree->computeLibpllLoglk() << endl;
        stats << "D " << jointTree->getDupRate() << endl;
        stats << "L " << jointTree->getLossRate() << endl;
        stats << "T " << jointTree->getTransferRate() << endl;
        stats << "hash " << jointTree->getUnrootedTreeHash() << endl;
        stats << " " << endl;
      }
      jointTree->save(allTreesFile, !firstRun);
      Scenario scenario;
      jointTree->inferMLScenario(scenario);
      Logger::info << endl;
      scenario.saveEventsCounts(eventCountsFile);
      scenario.saveTreeWithEvents(treeWithEventsFile);
    }
    firstRun = false;
  }  
  Logger::info << "Best tree: " + bestTreeFile << endl;
  Logger::info << "Best tree with events: " + treeWithEventsFile << endl;
  Logger::timed << "End of JointSearch execution" << endl;
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

