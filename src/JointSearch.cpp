#include "Arguments.hpp"
#include "ParallelContext.hpp"
#include <ale/tools/IO/IO.h>
#include <ale/tools/Utils.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/NNISearch.h>
#include <treeSearch/SPRSearch.h>
#include <Logger.hpp>
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
  double dupRate = 2;
  double lossRate = 1;
  
  Logger::init();
  Arguments::init(argc, argv);
  Arguments::checkInputs();
  Logger::initFileOutput(Arguments::output);
  Arguments::printCommand();
  Arguments::printSummary();
  auto jointTree = make_shared<JointTree>(Arguments::geneTree,
      Arguments::alignment,
      Arguments::speciesTree,
      dupRate,
      lossRate
      );
  printJointTreeInfo(jointTree);
  jointTree->optimizeParameters();
  Logger::timed << "Starting search..." << endl;
  if (Arguments::strategy == "SPR") {
    SPRSearch::applySPRSearch(*jointTree);
  } else if (Arguments::strategy == "NNI") {
    NNISearch::applyNNISearch(*jointTree);
  } else if (Arguments::strategy == "EVAL") {
  } else if (Arguments::strategy == "HYBRID") {
    NNISearch::applyNNISearch(*jointTree);
    jointTree->optimizeParameters();
    SPRSearch::applySPRSearch(*jointTree);
    jointTree->optimizeParameters();
    NNISearch::applyNNISearch(*jointTree);
  }
  Logger::timed << "End of search" << endl;
  jointTree->printLoglk();
  if (!ParallelContext::getRank()) {
    jointTree->save(Arguments::output + ".newick");
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

