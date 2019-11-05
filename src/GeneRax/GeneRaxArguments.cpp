
#include "GeneRaxArguments.hpp"
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <algorithm>
#include <vector>

GeneRaxArguments::GeneRaxArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  strategy(Strategy::SPR),
  reconciliationModelStr("UndatedDL"),
  reconciliationOpt(RecOpt::Grid),
  output("GeneRax"),
  perFamilyDTLRates(false),
  duplicates(1),
  initStrategies(3),
  rootedGeneTree(true),
  supportThreshold(-1.0),
  recRadius(0),
  perSpeciesDTLRates(false),
  useTransferFrequencies(false),
  userDTLRates(false),
  dupRate(1.0),
  lossRate(1.0),
  transferRate(0.0),
  optimizeGeneTrees(true),
  maxSPRRadius(5),
  recWeight(1.0), 
  seed(123),
  exec(iargv[0]),
  optimizeSpeciesTree(false),
  speciesFastRadius(5),
  speciesSlowRadius(0)
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
      strategy = ArgumentsHelper::strToStrategy(std::string(argv[++i]));
      if (strategy == Strategy::EVAL) {
        recRadius = maxSPRRadius = 0;
      }
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModelStr = std::string(argv[++i]);
    } else if (arg == "--rec-opt") {
      reconciliationOpt = ArgumentsHelper::strToRecOpt(std::string(argv[++i]));
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--per-family-rates") {
      perFamilyDTLRates = true;
    } else if (arg == "--init-strategies") {
      initStrategies = atoi(argv[++i]);
    } else if (arg == "--duplicates") {
      duplicates = atoi(argv[++i]);
    } else if (arg == "--unrooted-gene-tree") {
      rootedGeneTree = false;
    } else if (arg == "--support-threshold") {
      supportThreshold = static_cast<double>(atof(argv[++i]));
    } else if (arg == "--rec-radius") {
      recRadius = atoi(argv[++i]);
    } else if (arg == "--per-species-rates") {
      perSpeciesDTLRates = true;
    } else if (arg == "--use-transfer-frequencies") {
      useTransferFrequencies = true;
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
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--species-fast-radius") {
      speciesFastRadius = atoi(argv[++i]);
    } else if (arg == "--species-slow-radius") {
      speciesSlowRadius = atoi(argv[++i]);
    } else if (arg == "--optimize-species-tree") {
      optimizeSpeciesTree = true;
    } else if (arg == "--do-not-optimize-gene-trees") {
      optimizeGeneTrees = false;
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
  if (!speciesTree.size() && !optimizeSpeciesTree) {
    Logger::info << "[Error] You need to provide a species tree or to optimize it." << std::endl;
    ok = false;
  }
  if (userDTLRates && perSpeciesDTLRates) {
    Logger::info << "[Error] You cannot specify the rates when using per-species DTL rates" << std::endl;
    ok = false;
  }
  if (userDTLRates && (dupRate < 0.0 || lossRate < 0.0)) {
    Logger::info << "[Error] You specified at least one of the duplication and loss rates, but not both of them." << std::endl;
    ok = false;
  }
  if (perSpeciesDTLRates && perFamilyDTLRates) {
    Logger::info << "[Error] You cannot use per-family and per-species rates at the same time" << std::endl;
    ok = false;
  }
  if (!ArgumentsHelper::isValidRecModel(reconciliationModelStr)) {
    Logger::info << "[Error] Invalid reconciliation model string " << reconciliationModelStr << std::endl;
    ok = false;
  }
  if (!ok) {
    Logger::info << "Aborting." << std::endl;
    ParallelContext::abort(1);
  }
  if (speciesTree.size() && speciesTree != "random") {
    assertFileExists(speciesTree);
  }
}

void GeneRaxArguments::printHelp() {
  Logger::info << "-h, --help" << std::endl;
  Logger::info << "-f, --families <FAMILIES_INFORMATION>" << std::endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << std::endl;
  Logger::info << "--strategy <STRATEGY>  {EVAL, SPR}" << std::endl;
  Logger::info << "-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL, Auto}" << std::endl;
  Logger::info << "--rec-opt <reconciliationOpt>  {window, simplex}" << std::endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << std::endl;
  Logger::info << "--duplicates <DUPLICATES_NUMBER>" << std::endl;
  Logger::info << "--init-strategies <1 or 4>" << std::endl;
  Logger::info << "--unrooted-gene-tree" << std::endl;
  Logger::info << "--support-threshold <threshold>" << std::endl;
  Logger::info << "--per-family-rates" << std::endl;
  Logger::info << "--per-species-rates" << std::endl;
  Logger::info << "--dup-rate <duplication rate>" << std::endl;
  Logger::info << "--loss-rate <loss rate>" << std::endl;
  Logger::info << "--transfer-rate <transfer rate>" << std::endl;
  Logger::info << "--max-spr-radius <max SPR radius>" << std::endl;
  Logger::info << "--rec-weight <reconciliation likelihood weight>" << std::endl;
  Logger::info << "--seed <seed>" << std::endl;
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
  Logger::info << "Strategy: " << ArgumentsHelper::strategyToStr(strategy) << std::endl;
  Logger::info << "Reconciliation model: " << reconciliationModelStr << std::endl;
  Logger::info << "Reconciliation opt: " << ArgumentsHelper::recOptToStr(reconciliationOpt) << std::endl;
  Logger::info << "DTL rates: "; 
  if (perSpeciesDTLRates) {
    Logger::info << "per species rates" << std::endl;
  } else if (perFamilyDTLRates) {
    Logger::info << "per family rates" << std::endl;
  } else {
    Logger::info << "global rates" << std::endl;
  }
  Logger::info << "Prefix: " << output << std::endl;
  Logger::info << "Duplicates: " << duplicates << std::endl;
  if (duplicates > 1) {
    Logger::info << "Init strategies: " << initStrategies << std::endl;
  }
  Logger::info << "Unrooted gene tree: " << boolStr[!rootedGeneTree] << std::endl;
  Logger::info << "Reconciliation radius: " << recRadius << std::endl;
#ifdef WITH_MPI
  Logger::info << "MPI Ranks: " << ParallelContext::getSize() << std::endl;
#else
  Logger::info << "You are running GeneRax without MPI (no parallelization)" << std::endl;
#endif
  Logger::info << "Max SPR radius: " << maxSPRRadius << std::endl;
  Logger::info << "Support threshold: " << supportThreshold << std::endl;
  Logger::info << "Reconciliation likelihood weight: " << recWeight << std::endl;
  Logger::info << "Random seed: " << seed << std::endl;
  Logger::info << std::endl;
}
