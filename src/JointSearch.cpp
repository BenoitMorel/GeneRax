#include <ale/tools/IO/IO.h>
#include <ale/tools/Utils.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/NNISearch.h>
#include <treeSearch/SPRSearch.h>
#include <chrono>

using namespace std;

struct Arguments {
  Arguments(int argc, char * argv[]) {
    if (argc != 7) {
      cerr << "argc = " << argc << endl;
      for (int i = 0; i < argc; ++i) {
        cerr << argv[i] << " ";
      }
      cerr << endl;
      cerr << "Syntax error. Usage:" << endl;
      cerr << "./jointTreeSearch gene_tree alignment_file species_tree strategy threads output_file" << endl;
      cerr << "Strategy can be NNI, SPR or HYBRID" << endl;
      exit(1);
    }
    int i = 1;
    geneTree = string(argv[i++]);
    alignment = string(argv[i++]);
    speciesTree = string(argv[i++]);
    strategy = string(argv[i++]);
    threads = atoi(argv[i++]);
    output = string(argv[i++]);

  }
  
  string geneTree;
  string alignment;
  string speciesTree;
  string strategy;
  int threads;
  string output;
};


int main(int argc, char * argv[]) {
  auto start = chrono::high_resolution_clock::now();
  double dupCost = 2;
  double lossCost = 1;
  Arguments arg(argc, argv);
  shared_ptr<AbstractJointTree> jointTree = make_shared<ParallelJointTree>(arg.geneTree,
      arg.alignment,
      arg.speciesTree,
      dupCost,
      lossCost,
      arg.threads
      );
  jointTree->optimizeParameters();
  if (arg.strategy == "SPR") {
    cout << "Starting SPR search" << endl;
    SPRSearch::applySPRSearch(*jointTree);
  } else if (arg.strategy == "NNI") {
    cout << "Starting NNI search" << endl;
    NNISearch::applyNNISearch(*jointTree);
  } else if (arg.strategy == "HYBRID") {
    cout << "Starting NNI search" << endl;
    NNISearch::applyNNISearch(*jointTree);
    cout << "Starting SPR search" << endl;
    SPRSearch::applySPRSearch(*jointTree);
    cout << "Starting NNI search" << endl;
    NNISearch::applyNNISearch(*jointTree);
  }
  cout << "END OF THE SEARCH" << endl;

  jointTree->getThreadInstance().printLoglk();
  jointTree->getThreadInstance().save(arg.output + ".newick");

  auto finish = chrono::high_resolution_clock::now();
  
  ofstream os(arg.output + ".stats");
  jointTree->getThreadInstance().printLoglk(true, true, true, os);

  chrono::duration<double> elapsed = finish - start;
  cout<<"Elapsed time is :  "<< chrono::duration_cast<chrono::seconds>(elapsed).count() << "s" <<endl;  
  return 0;
}

