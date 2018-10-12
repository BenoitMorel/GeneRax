#include "Arguments.hpp"
#include <Logger.hpp>

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
      Logger::error << "Unrecognized argument " << arg << endl;
      Logger::error << "Aborting" << endl;
      exit(1);
    }
  }
}

void Arguments::checkInputs() {
  bool ok = true;
  if (!geneTree.size()) {
    Logger::error << "You need to provide a gene tree." << endl;
    ok = false;
  }
  if (!alignment.size()) {
    Logger::error << "You need to provide an alignment." << endl;
    ok = false;
  }
  if (!speciesTree.size()) {
    Logger::error << "You need to provide a species tree." << endl;
    ok = false;
  }
  if (!ok) {
    Logger::error << "Aborting." << endl;
    exit(1);
  }
}

void Arguments::printHelp() {
  Logger::error << "-h, --help" << endl;
  Logger::error << "-g, --gene-tree <GENE TREE>" << endl;
  Logger::error << "-a, --alignment <ALIGNMENT>" << endl;
  Logger::error << "-s, --species-tree <SPECIES TREE>" << endl;
  Logger::error << "--strategy <STRATEGY>" << endl;
  Logger::error << "-t, --threads <THREADS NUMBER>" << endl;
  Logger::error << "-p, --prefix <OUTPUT PREFIX>" << endl;
  Logger::error << "--check" << endl;
  Logger::error << "--verbose" << endl;
  Logger::error << endl;

}

void Arguments::printCommand() {
  Logger::info << "JointSearch was called as follow:" << endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << endl << endl;
}

void Arguments::printSummary() {
  string boolStr[2] = {string("OFF"), string("ON")};
  Logger::info << "Parameters summary: " << endl;
  Logger::info << "Gene tree: " << geneTree << endl;
  Logger::info << "Alignment: " << alignment << endl; 
  Logger::info << "Species tree: " << speciesTree << endl;
  Logger::info << "Strategy: " << strategy << endl;
  Logger::info << "Prefix: " << output << endl;
  Logger::info << "Cores: " << threads << endl;
  Logger::info << "Check mode: " << boolStr[check] << endl;
  Logger::info << "Verbose mode: " << boolStr[verbose] << endl;
  Logger::info << endl;
}
