
#include "GeneRaxArguments.hpp"
#include <IO/Logger.hpp>
#include <ParallelContext.hpp>
#include <algorithm>
#include <vector>

GeneRaxArguments::GeneRaxArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  strategy(SPR),
  reconciliationModelStr("UndatedDL"),
  reconciliationOpt(Simplex),
  output("GeneRax"),
  rootedGeneTree(true),
  perSpeciesDTLRates(false),
  userDTLRates(false),
  dupRate(1.0),
  lossRate(1.0),
  transferRate(0.0),
  maxSPRRadius(5),
  recWeight(1.0)
{
  if (argc == 1) {
    printHelp();
    ParallelContext::abort(0);
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      ParallelContext::abort(0);
    } else if (arg == "-f" || arg == "--families") {
      families = std::string(argv[++i]);
    } else if (arg == "-s" || arg == "--species-tree") {
      speciesTree = std::string(argv[++i]);
    } else if (arg == "--strategy") {
      strategy = Arguments::strToStrategy(std::string(argv[++i]));
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModelStr = std::string(argv[++i]);
      //reconciliationModel = Arguments::strToRecModel(std::string(argv[++i]));
    } else if (arg == "--rec-opt") {
      reconciliationOpt = Arguments::strToRecOpt(std::string(argv[++i]));
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--unrooted-gene-tree") {
      rootedGeneTree = false;
    } else if (arg == "--per-species-rates") {
      perSpeciesDTLRates = true;
    } else if (arg == "--dup-rate") {
      dupRate = atof(argv[++i]);
      userDTLRates = true;
    } else if (arg == "--loss-rate") {
      lossRate = atof(argv[++i]);
      userDTLRates = true;
    } else if (arg == "--transfer-rate") {
      transferRate = atof(argv[++i]);
      userDTLRates = true;
    } else if (arg == "--max-spr-radius") {
      maxSPRRadius = atoi(argv[++i]);
    } else if (arg == "--rec-weight") {
      recWeight = atof(argv[++i]);
    } else {
      Logger::error << "Unrecognized argument " << arg << std::endl;
      Logger::error << "Aborting" << std::endl;
      ParallelContext::abort(1);
    }
  }
  execPath = std::string(argv[0]);
  checkInputs();
}

void assertFileExists(const std::string &file) 
{
  std::ifstream f(file);
  if (!f) {
    Logger::error << "File " << file << " does not exist. Aborting." << std::endl;
    ParallelContext::abort(1);
  }
}

bool isIn(const std::string &elem, const std::vector<std::string> &v) {
  return find(v.begin(), v.end(), elem) != v.end();
}

void GeneRaxArguments::checkInputs() {
  bool ok = true;
  if (!speciesTree.size()) {
    Logger::error << "You need to provide a species tree." << std::endl;
    ok = false;
  }
  if (userDTLRates && perSpeciesDTLRates) {
    Logger::error << "You cannot specify the rates when using per-species DTL rates" << std::endl;
    ok = false;
  }
  if (userDTLRates && (dupRate < 0.0 || lossRate < 0.0)) {
    Logger::error << "You specified at least one of the duplication and loss rates, but not both of them." << std::endl;
    ok = false;
  }
  if (userDTLRates && reconciliationModelStr == "AutoDetect") {
    Logger::error << "You cannot set the DTL rates in AutoDetect mode" << std::endl;
    ok = false;
  }
  if (!ok) {
    Logger::error << "Aborting." << std::endl;
    ParallelContext::abort(1);
  }
  
  assertFileExists(speciesTree);
}

void GeneRaxArguments::printHelp() {
  Logger::info << "-h, --help" << std::endl;
  Logger::info << "-f, --families <FAMILIES_INFORMATION>" << std::endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << std::endl;
  Logger::info << "--strategy <STRATEGY>  {EVAL, SPR}" << std::endl;
  Logger::info << "-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL, Auto}" << std::endl;
  Logger::info << "--rec-opt <reconciliationOpt>  {window, simplex}" << std::endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << std::endl;
  Logger::info << "--unrooted-gene-tree" << std::endl;
  Logger::info << "--per-species-rates" << std::endl;
  Logger::info << "--dup-rate <duplication rate>" << std::endl;
  Logger::info << "--loss-rate <loss rate>" << std::endl;
  Logger::info << "--transfer-rate <transfer rate>" << std::endl;
  Logger::info << "--max-spr-radius <max SPR radius>" << std::endl;
  Logger::info << "--rec-weight <reconciliation likelihood weight>" << std::endl;
  Logger::info << "Please find more information on the GeneRax github wiki" << std::endl;
  Logger::info << std::endl;

}

void GeneRaxArguments::printCommand() {
  Logger::info << "GeneRax was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

void GeneRaxArguments::printSummary() {
  std::string boolStr[2] = {std::string("OFF"), std::string("ON")};
  Logger::info << "Parameters summary: " << std::endl;
  Logger::info << "Families information: " << families << std::endl;
  Logger::info << "Species tree: " << speciesTree << std::endl;
  Logger::info << "Strategy: " << Arguments::strategyToStr(strategy) << std::endl;
  Logger::info << "Reconciliation model: " << reconciliationModelStr << std::endl;
  Logger::info << "Reconciliation opt: " << Arguments::recOptToStr(reconciliationOpt) << std::endl;
  Logger::info << "DTL rates: " << (perSpeciesDTLRates ? "per-species" : "global") << std::endl;
  Logger::info << "Prefix: " << output << std::endl;
  Logger::info << "Unrooted gene tree: " << boolStr[!rootedGeneTree] << std::endl;
  Logger::info << "MPI Ranks: " << ParallelContext::getSize() << std::endl;
  Logger::info << "Max SPR radius: " << maxSPRRadius << std::endl;
  Logger::info << "Reconciliation likelihood weight: " << recWeight << std::endl;
  Logger::info << std::endl;
}
