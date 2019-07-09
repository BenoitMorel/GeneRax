
#include "SpeciesRaxArguments.hpp"
#include <IO/Logger.hpp>
#include "Arguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <algorithm>
#include <vector>

SpeciesRaxArguments::SpeciesRaxArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  reconciliationModel(UndatedDL),
  reconciliationOpt(Simplex),
  output("SpeciesRax"),
  rootedGeneTree(true),
  userDTLRates(false),
  dupRate(1.0),
  lossRate(1.0),
  transferRate(0.0)
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
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModel = Arguments::strToRecModel(std::string(argv[++i]));
    } else if (arg == "--rec-opt") {
      reconciliationOpt = Arguments::strToRecOpt(std::string(argv[++i]));
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
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
    } else {
      Logger::error << "Unrecognized argument " << arg << std::endl;
      Logger::error << "Aborting" << std::endl;
      ParallelContext::abort(1);
    }
  }
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

void SpeciesRaxArguments::checkInputs() {
  bool ok = true;
  if (!speciesTree.size()) {
    Logger::error << "You need to provide a species tree." << std::endl;
    ok = false;
  }
  if (userDTLRates && (dupRate < 0.0 || lossRate < 0.0)) {
    Logger::error << "You specified at least one of the duplication and loss rates, but not both of them." << std::endl;
    ok = false;
  }
  if (!ok) {
    Logger::error << "Aborting." << std::endl;
    ParallelContext::abort(1);
  }
  
  assertFileExists(speciesTree);
}

void SpeciesRaxArguments::printHelp() {
  Logger::info << "-h, --help" << std::endl;
  Logger::info << "-f, --families <FAMILIES_INFORMATION>" << std::endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << std::endl;
  Logger::info << "-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL}" << std::endl;
  Logger::info << "--rec-opt <reconciliationOpt>  {window, simplex}" << std::endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << std::endl;
  Logger::info << "--unrooted-gene-tree" << std::endl;
  Logger::info << "--dupRate <duplication rate>" << std::endl;
  Logger::info << "--lossRate <loss rate>" << std::endl;
  Logger::info << "--transferRate <transfer rate>" << std::endl;
  Logger::info << std::endl;

}

void SpeciesRaxArguments::printCommand() {
  Logger::info << "SpeciesRax was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

void SpeciesRaxArguments::printSummary() {
  std::string boolStr[2] = {std::string("OFF"), std::string("ON")};
  Logger::info << "Parameters summary: " << std::endl;
  Logger::info << "Families information: " << families << std::endl;
  Logger::info << "Species tree: " << speciesTree << std::endl;
  Logger::info << "Reconciliation model: " << Arguments::recModelToStr(reconciliationModel) << std::endl;
  Logger::info << "Reconciliation opt: " << Arguments::recOptToStr(reconciliationOpt) << std::endl;
  Logger::info << "Prefix: " << output << std::endl;
  Logger::info << "Unrooted gene tree: " << boolStr[!rootedGeneTree] << std::endl;
  Logger::info << "MPI Ranks: " << ParallelContext::getSize() << std::endl;
  Logger::info << std::endl;
}
