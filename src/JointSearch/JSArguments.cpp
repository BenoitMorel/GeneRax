#include "JSArguments.hpp"
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <algorithm>
#include <vector>
#include <IO/ArgumentsHelper.hpp>

JSArguments::JSArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  strategy(GeneSearchStrategy::EVAL),
  reconciliationModel(RecModel::UndatedDL),
  reconciliationOpt(RecOpt::Simplex),
  libpllModel("GTR"),
  output("jointSearch"),
  check(false),
  rootedGeneTree(true),
  userDTLRates(false),
  dupRate(-1.0),
  lossRate(-1.0),
  transferRate(-1.0)
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
    } else if (arg == "-g" || arg == "--gene-tree") {
      geneTree = std::string(argv[++i]);
    } else if (arg == "-a" || arg == "--alignment") {
      alignment = std::string(argv[++i]);
    } else if (arg == "-s" || arg == "--species-tree") {
      speciesTree = std::string(argv[++i]);
    } else if (arg == "-m" || arg == "--map") {
      geneSpeciesMap = std::string(argv[++i]);
    } else if (arg == "--strategy") {
      strategy = ArgumentsHelper::strToStrategy(std::string(argv[++i]));
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModel = ArgumentsHelper::strToRecModel(std::string(argv[++i]));
    } else if (arg == "--rec-opt") {
      reconciliationOpt = ArgumentsHelper::strToRecOpt(std::string(argv[++i]));
    } else if (arg == "--libpll-model") {
      libpllModel = std::string(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--check") {
      check = true;
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

static void assertFileExists(const std::string &file) 
{
  std::ifstream f(file);
  if (!f) {
    Logger::error << "File " << file << " does not exist. Aborting." << std::endl;
    ParallelContext::abort(1);
  }
}

void JSArguments::checkInputs() {
  bool ok = true;
  if (!alignment.size()) {
    Logger::error << "You need to provide an alignment." << std::endl;
    ok = false;
  }
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
  
  if (geneTree.size() && geneTree != "__random__") {
    assertFileExists(geneTree);
  }
  assertFileExists(speciesTree);
  if (geneSpeciesMap.size()) {
    assertFileExists(geneSpeciesMap);
  }
  assertFileExists(alignment);
}

void JSArguments::printHelp() {
  Logger::info << "-h, --help" << std::endl;
  Logger::info << "-g, --gene-tree <GENE TREE>" << std::endl;
  Logger::info << "-a, --alignment <ALIGNMENT>" << std::endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << std::endl;
  Logger::info << "-m, --map <GENE_SPECIES_MAPPING>" << std::endl;
  Logger::info << "--strategy <STRATEGY>  {EVAL, SPR}" << std::endl;
  Logger::info << "-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL}" << std::endl;
  Logger::info << "--rec-opt <reconciliationOpt>  {grid, simplex}" << std::endl;
  Logger::info << "--libpll-model <libpllModel>  {GTR, LG, DAYHOFF etc.}" << std::endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << std::endl;
  Logger::info << "--check" << std::endl;
  Logger::info << "--unrooted-gene-tree" << std::endl;
  Logger::info << "--dupRate <duplication rate>" << std::endl;
  Logger::info << "--lossRate <loss rate>" << std::endl;
  Logger::info << "--transferRate <transfer rate>" << std::endl;
  Logger::info << std::endl;

}

void JSArguments::printCommand() {
  Logger::info << "JointSearch was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

void JSArguments::printSummary() {
  std::string boolStr[2] = {std::string("OFF"), std::string("ON")};
  Logger::info << "Parameters summary: " << std::endl;
  Logger::info << "Gene tree: " << geneTree << std::endl;
  Logger::info << "Alignment: " << alignment << std::endl; 
  Logger::info << "Species tree: " << speciesTree << std::endl;
  Logger::info << "Gene species map: " << geneSpeciesMap << std::endl;
  Logger::info << "Strategy: " << ArgumentsHelper::strategyToStr(strategy) << std::endl;
  Logger::info << "Reconciliation model: " << ArgumentsHelper::recModelToStr(reconciliationModel) << std::endl;
  Logger::info << "Reconciliation opt: " << ArgumentsHelper::recOptToStr(reconciliationOpt) << std::endl;
  Logger::info << "Libpll model: " << libpllModel << std::endl;
  Logger::info << "Prefix: " << output << std::endl;
  Logger::info << "Check mode: " << boolStr[check] << std::endl;
  Logger::info << "Unrooted gene tree: " << boolStr[!rootedGeneTree] << std::endl;
  Logger::info << "MPI Ranks: " << ParallelContext::getSize() << std::endl;
  Logger::info << std::endl;
}
