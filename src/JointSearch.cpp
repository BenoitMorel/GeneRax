#include "Arguments.hpp"
#include "ParallelContext.hpp"
#include <treeSearch/JointTree.h>
#include <treeSearch/SPRSearch.h>
#include <Logger.hpp>
#include <algorithm>
#include <limits>
using namespace std;


void printJointTreeInfo(shared_ptr<JointTree> tree) 
{
  auto treeInfo = tree->getTreeInfo();
  int speciesLeaves = tree->getSpeciesTree()->tip_count;
  int geneLeaves = treeInfo->tip_count;;
  int sites = treeInfo->partitions[0]->sites;


  Logger::info << "Species leaves: " << speciesLeaves << endl;
  Logger::info << "Gene leaves: " << geneLeaves << endl;
  Logger::info << "Sites: " << sites << endl;
  Logger::info << endl;
}



int internal_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  double dupRate = 0.2;
  double lossRate = 0.1;
  
  Logger::init();
  Arguments::init(argc, argv);
  Arguments::checkInputs();
  Logger::initFileOutput(Arguments::output);
  Arguments::printCommand();
  Arguments::printSummary();
  
  vector<string> geneTreeStrings;
  string geneTreeString;
  ifstream treeStream(Arguments::geneTree);
  while(getline(treeStream, geneTreeString)) {
    geneTreeString.erase(remove(geneTreeString.begin(), geneTreeString.end(), '\n'), geneTreeString.end());
    if (geneTreeString.empty()) {
      continue;
    }
    geneTreeStrings.push_back(geneTreeString);
  }
  bool firstRun  = true; 
  double bestLL = numeric_limits<double>::lowest();
  for (auto &geneTreeString: geneTreeStrings) {
    auto jointTree = make_shared<JointTree>(geneTreeString,
        Arguments::alignment,
        Arguments::speciesTree,
        Arguments::geneSpeciesMap,
        Arguments::reconciliationModel,
        dupRate,
        lossRate
        );
    
    
    printJointTreeInfo(jointTree);
    jointTree->optimizeParameters();
    Logger::timed << "Starting search..." << endl;
    if (Arguments::strategy == "SPR") {
      SPRSearch::applySPRSearch(*jointTree);
    } else if (Arguments::strategy == "EVAL") {
    }
    Logger::timed << "End of search" << endl;
    jointTree->printLoglk();
    Logger::info << "Final tree hash: " << jointTree->getUnrootedTreeHash() << endl;
    if (!ParallelContext::getRank()) {
      double ll = jointTree->computeJointLoglk();
      assert(!isnan(ll));
      if (ll >= bestLL) {
        bestLL = ll;
        jointTree->save(Arguments::output + ".newick", false);
      }
      jointTree->save(Arguments::output + "_all" + ".newick", !firstRun);
    }
    firstRun = false;
  }  
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

