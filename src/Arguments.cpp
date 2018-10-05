#include "Arguments.hpp"
#include <iostream>

int Arguments::argc = 0;
char ** Arguments::argv = 0;
string Arguments::geneTree;
string Arguments::alignment;
string Arguments::speciesTree;
string Arguments::strategy("NNI");
int Arguments::threads = 1;
string Arguments::output("jointSearch");
bool Arguments::check = false;
bool Arguments::verbose = false;
double Arguments::aleWeight = 1.0;
bool Arguments::costsEstimation = false;

void Arguments::init(int argc, char * argv[])
{
  Arguments::argc = argc;
  Arguments::argv = argv;
  if (argc == 1) {
    printHelp();
    exit(0);
  }
  for (int i = 1; i < argc; ++i) {
    string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      exit(0);
    } else if (arg == "-g" || arg == "--gene-tree") {
      geneTree = string(argv[++i]);
    } else if (arg == "-a" || arg == "--alignment") {
      alignment = string(argv[++i]);
    } else if (arg == "-s" || arg == "--species-tree") {
      speciesTree = string(argv[++i]);
    } else if (arg == "--strategy") {
      strategy = string(argv[++i]);
    } else if (arg == "-t" || arg == "--threads") {
      threads = atoi(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      output = string(argv[++i]);
    } else if (arg == "-c" || arg == "--cost-estimation") {
      costsEstimation = true;
    } else if (arg == "--check") {
      check = true;
    } else if (arg == "--verbose") {
      verbose = true;
    } else if (arg == "--ale-weight") {
      aleWeight = atof(argv[++i]);
    } else {
      cerr << "Unrecognized argument " << arg << endl;
      cerr << "Aborting" << endl;
      exit(1);
    }
  }
}

void Arguments::checkInputs() {
  bool ok = true;
  if (!geneTree.size()) {
    cerr << "You need to provide a gene tree." << endl;
    ok = false;
  }
  if (!alignment.size()) {
    cerr << "You need to provide an alignment." << endl;
    ok = false;
  }
  if (!speciesTree.size()) {
    cerr << "You need to provide a species tree." << endl;
    ok = false;
  }
  if (!ok) {
    cerr << "Aborting." << endl;
    exit(1);
  }
}

void Arguments::printHelp() {
  cerr << "-h, --help" << endl;
  cerr << "-g, --gene-tree <GENE TREE>" << endl;
  cerr << "-a, --alignment <ALIGNMENT>" << endl;
  cerr << "-s, --species-tree <SPECIES TREE>" << endl;
  cerr << "--strategy <STRATEGY>" << endl;
  cerr << "-t, --threads <THREADS NUMBER>" << endl;
  cerr << "-p, --prefix <OUTPUT PREFIX>" << endl;
  cerr << "--check" << endl;
  cerr << "--verbose" << endl;

}

void Arguments::printCommand() {
  cout << "JointSearch was called as follow:" << endl;
  for (int i = 0; i < argc; ++i) {
    cout << argv[i] << " ";
  }
  cout << endl;
}

void Arguments::printSummary() {
  string boolStr[2] = {string("OFF"), string("ON")};
  cout << "Parameters summary: " << endl;
  cout << "Gene tree: " << geneTree << endl;
  cout << "Alignment: " << alignment << endl; 
  cout << "Species tree: " << speciesTree << endl;
  cout << "Strategy: " << strategy << endl;
  cout << "Prefix: " << output << endl;
  cout << "Threads: " << threads << endl;
  cout << "Check mode: " << boolStr[check] << endl;
  cout << "Verbose mode: " << boolStr[verbose] << endl;
}
