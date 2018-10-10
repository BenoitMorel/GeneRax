#include "Arguments.hpp"
#include "ParallelContext.hpp"
#include <ale/tools/IO/IO.h>
#include <ale/tools/Utils.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/NNISearch.h>
#include <treeSearch/SPRSearch.h>
#include <chrono>

using namespace std;


int main(int argc, char * argv[]) {
  ParallelContext::init();
  auto start = chrono::high_resolution_clock::now();
  double dupRate = 2;
  double lossRate = 1;
  Arguments::init(argc, argv);
  Arguments::printCommand();
  cout << endl;
  Arguments::printSummary();
  cout << endl;
  auto jointTree = make_shared<JointTree>(Arguments::geneTree,
      Arguments::alignment,
      Arguments::speciesTree,
      dupRate,
      lossRate
      );
  jointTree->optimizeParameters();
  if (Arguments::costsEstimation) {
    cout << "Starting cost estimation..." << endl;
    jointTree->optimizeDTRates();
  }
  if (Arguments::strategy == "SPR") {
    cout << "Starting SPR search" << endl;
    SPRSearch::applySPRSearch(*jointTree);
  } else if (Arguments::strategy == "NNI") {
    cout << "Starting NNI search" << endl;
    NNISearch::applyNNISearch(*jointTree);
  } else if (Arguments::strategy == "EVAL") {
  } else if (Arguments::strategy == "HYBRID") {
    cout << "Starting NNI search" << endl;
    NNISearch::applyNNISearch(*jointTree);
    cout << "Starting SPR search" << endl;
    SPRSearch::applySPRSearch(*jointTree);
    cout << "Starting NNI search" << endl;
    NNISearch::applyNNISearch(*jointTree);
  }

  jointTree->printLoglk();
  jointTree->save(Arguments::output + ".newick");

  auto finish = chrono::high_resolution_clock::now();
  
  ofstream os(Arguments::output + ".stats");
  jointTree->printLoglk(true, true, true, os);

  chrono::duration<double> elapsed = finish - start;
  cout<<"Elapsed time is :  "<< chrono::duration_cast<chrono::milliseconds>(elapsed).count() / 1000.0 << "s" <<endl;  
  ParallelContext::finalize();
  return 0;
}

