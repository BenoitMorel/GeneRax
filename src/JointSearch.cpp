#include "Arguments.hpp"
#include "ParallelContext.hpp"
#include <ale/tools/IO/IO.h>
#include <ale/tools/Utils.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/NNISearch.h>
#include <treeSearch/SPRSearch.h>
#include <Logger.hpp>
using namespace std;


int internal_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  Logger::init();
  double dupRate = 2;
  double lossRate = 1;
  
  Arguments::init(argc, argv);
  Arguments::printCommand();
  Arguments::printSummary();
  auto jointTree = make_shared<JointTree>(Arguments::geneTree,
      Arguments::alignment,
      Arguments::speciesTree,
      dupRate,
      lossRate
      );
  jointTree->optimizeParameters();
  if (Arguments::costsEstimation) {
    Logger::timed << "Starting cost estimation..." << endl;
    jointTree->optimizeDTRates();
  }
  Logger::timed << "Starting search..." << endl;
  if (Arguments::strategy == "SPR") {
    SPRSearch::applySPRSearch(*jointTree);
  } else if (Arguments::strategy == "NNI") {
    NNISearch::applyNNISearch(*jointTree);
  } else if (Arguments::strategy == "EVAL") {
  } else if (Arguments::strategy == "HYBRID") {
    NNISearch::applyNNISearch(*jointTree);
    SPRSearch::applySPRSearch(*jointTree);
    NNISearch::applyNNISearch(*jointTree);
  }
  Logger::timed << "End of search" << endl;
  jointTree->printLoglk();
  jointTree->save(Arguments::output + ".newick");

  ofstream os(Arguments::output + ".stats");
  jointTree->printLoglk(true, true, true, os);

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

