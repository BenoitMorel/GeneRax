
#include "GeneRaxArguments.hpp"
#include <IO/Logger.hpp>
#include "Arguments.hpp"
#include <ParallelContext.hpp>
#include <algorithm>
#include <vector>

GeneRaxArguments::GeneRaxArguments(int argc, char * argv[]):
  argc(argc),
  argv(argv),
  strategy(EVAL),
  reconciliationModel(UndatedDL),
  reconciliationOpt(Simplex),
  output("GeneRax"),
  rootedGeneTree(true),
  userDTLRates(false),
  dupRate(1.0),
  lossRate(1.0),
  transferRate(0.0),
  autodetectDTLModel(false),
  maxSPRRadius(5)
{
  if (argc == 1) {
    printHelp();
    ParallelContext::abort(0);
  }
  for (int i = 1; i < argc; ++i) {
    string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      ParallelContext::abort(0);
    } else if (arg == "-f" || arg == "--families") {
      families = string(argv[++i]);
    } else if (arg == "-s" || arg == "--species-tree") {
      speciesTree = string(argv[++i]);
    } else if (arg == "--strategy") {
      strategy = Arguments::strToStrategy(string(argv[++i]));
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModel = Arguments::strToRecModel(string(argv[++i]));
    } else if (arg == "--rec-opt") {
      reconciliationOpt = Arguments::strToRecOpt(string(argv[++i]));
    } else if (arg == "-p" || arg == "--prefix") {
      output = string(argv[++i]);
    } else if (arg == "--unrooted-gene-tree") {
      rootedGeneTree = false;
    } else if (arg == "--dupRate") {
      dupRate = atof(argv[++i]);
      userDTLRates = true;
    } else if (arg == "--lossRate") {
      lossRate = atof(argv[++i]);
      userDTLRates = true;
    } else if (arg == "--transferRate") {
      transferRate = atof(argv[++i]);
      userDTLRates = true;
    } else if (arg == "--autodetectDTLRates") {
      autodetectDTLModel = true;
    } else if (arg == "--max-spr-radius") {
      maxSPRRadius = atoi(argv[++i]);
    } else {
      Logger::error << "Unrecognized argument " << arg << endl;
      Logger::error << "Aborting" << endl;
      ParallelContext::abort(1);
    }
  }
  checkInputs();
}

void assertFileExists(const string &file) 
{
  ifstream f(file);
  if (!f) {
    Logger::error << "File " << file << " does not exist. Aborting." << endl;
    ParallelContext::abort(1);
  }
}

bool isIn(const string &elem, const vector<string> &v) {
  return find(v.begin(), v.end(), elem) != v.end();
}

void GeneRaxArguments::checkInputs() {
  bool ok = true;
  if (!speciesTree.size()) {
    Logger::error << "You need to provide a species tree." << endl;
    ok = false;
  }
  if (userDTLRates && (dupRate < 0.0 || lossRate < 0.0)) {
    Logger::error << "You specified at least one of the duplication and loss rates, but not both of them." << endl;
    ok = false;
  }
  if (!ok) {
    Logger::error << "Aborting." << endl;
    ParallelContext::abort(1);
  }
  
  assertFileExists(speciesTree);
}

void GeneRaxArguments::printHelp() {
  Logger::info << "-h, --help" << endl;
  Logger::info << "-f, --families <FAMILIES_INFORMATION>" << endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << endl;
  Logger::info << "--strategy <STRATEGY>  {EVAL, SPR}" << endl;
  Logger::info << "-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL, DatedDL}" << endl;
  Logger::info << "--rec-opt <reconciliationOpt>  {window, simplex}" << endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << endl;
  Logger::info << "--unrooted-gene-tree" << endl;
  Logger::info << "--dupRate <duplication rate>" << endl;
  Logger::info << "--lossRate <loss rate>" << endl;
  Logger::info << "--transferRate <transfer rate>" << endl;
  Logger::info << "--max-spr-radius <max SPR radius>" << endl;
  Logger::info << endl;

}

void GeneRaxArguments::printCommand() {
  Logger::info << "GeneRax was called as follow:" << endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << endl << endl;
}

void GeneRaxArguments::printSummary() {
  string boolStr[2] = {string("OFF"), string("ON")};
  Logger::info << "Parameters summary: " << endl;
  Logger::info << "Families information: " << families << endl;
  Logger::info << "Species tree: " << speciesTree << endl;
  Logger::info << "Strategy: " << Arguments::strategyToStr(strategy) << endl;
  Logger::info << "Reconciliation model: " << Arguments::recModelToStr(reconciliationModel) << endl;
  Logger::info << "Reconciliation opt: " << Arguments::recOptToStr(reconciliationOpt) << endl;
  Logger::info << "Prefix: " << output << endl;
  Logger::info << "Unrooted gene tree: " << boolStr[!rootedGeneTree] << endl;
  Logger::info << "MPI Ranks: " << ParallelContext::getSize() << endl;
  Logger::info << "Max SPR radius: " << maxSPRRadius << endl;
  Logger::info << endl;
}
