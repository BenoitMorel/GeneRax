#include "GecoArguments.hpp"
#include <IO/Logger.hpp>
#include <ParallelContext.hpp>
#include <algorithm>
#include <vector>

GecoArguments::GecoArguments(int argc, char * argv[]):
  argc(argc),
  argv(argv),
  reconciliationModel("UndatedDL"),
  reconciliationOpt("simplex"),
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
    } else if (arg == "-m" || arg == "--map") {
      geneSpeciesMap = string(argv[++i]);
    } else if (arg == "--strategy") {
      strategy = string(argv[++i]);
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModel = string(argv[++i]);
    } else if (arg == "--rec-opt") {
      reconciliationOpt = string(argv[++i]);
    } else if (arg == "--libpll-model") {
      libpllModel = string(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      output = string(argv[++i]);
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
      Logger::error << "Unrecognized argument " << arg << endl;
      Logger::error << "Aborting" << endl;
      exit(1);
    }
  }
  checkInputs();
}

void assertFileExists(const string &file) 
{
  ifstream f(file);
  if (!f) {
    Logger::error << "File " << file << " does not exist. Aborting." << endl;
    exit(1);
  }
}

bool isIn(const string &elem, const vector<string> &v) {
  return find(v.begin(), v.end(), elem) != v.end();
}

void GecoArguments::checkInputs() {
  bool ok = true;
  if (!alignment.size()) {
    Logger::error << "You need to provide an alignment." << endl;
    ok = false;
  }
  if (!speciesTree.size()) {
    Logger::error << "You need to provide a species tree." << endl;
    ok = false;
  }
  if (!geneSpeciesMap.size()) {
    Logger::error << "You need to provide a gene species map file." << endl;
    ok = false;
  }
  if (userDTLRates && (dupRate < 0.0 || lossRate < 0.0)) {
    Logger::error << "You specified at least one of the duplication and loss rates, but not both of them." << endl;
    ok = false;
  }
  vector<string> possibleReconciliationModels;
  possibleReconciliationModels.push_back("UndatedDL");
  possibleReconciliationModels.push_back("UndatedDTL");
  possibleReconciliationModels.push_back("DatedDL");
  if (!isIn(reconciliationModel, possibleReconciliationModels)) {
    Logger::error << "Invalid reconciliation model " << reconciliationModel << endl;
    ok = false;
  }
  vector <string> possibleRecOpts;
  possibleRecOpts.push_back("window");
  possibleRecOpts.push_back("simplex");
  if (!isIn(reconciliationOpt, possibleRecOpts)) {
    Logger::error << "Invalid reconciliation opt " << reconciliationOpt << endl;
    ok = false;
  }
  if (!ok) {
    Logger::error << "Aborting." << endl;
    exit(1);
  }
  
  if (geneTree.size() && geneTree != "__random__") {
    assertFileExists(geneTree);
  }
  assertFileExists(speciesTree);
  assertFileExists(geneSpeciesMap);
  assertFileExists(alignment);
}

void GecoArguments::printHelp() {
  Logger::info << "-h, --help" << endl;
  Logger::info << "-g, --gene-tree <GENE TREE>" << endl;
  Logger::info << "-a, --alignment <ALIGNMENT>" << endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << endl;
  Logger::info << "-m, --map <GENE_SPECIES_MAPPING>" << endl;
  Logger::info << "--strategy <STRATEGY>  {EVAL, SPR}" << endl;
  Logger::info << "-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL, DatedDL}" << endl;
  Logger::info << "--rec-opt <reconciliationOpt>  {window, simplex}" << endl;
  Logger::info << "--libpll-model <libpllModel>  {GTR, LG, DAYHOFF etc.}" << endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << endl;
  Logger::info << "--check" << endl;
  Logger::info << "--unrooted-gene-tree" << endl;
  Logger::info << "--dupRate <duplication rate>" << endl;
  Logger::info << "--lossRate <loss rate>" << endl;
  Logger::info << "--transferRate <transfer rate>" << endl;
  Logger::info << endl;

}

void GecoArguments::printCommand() {
  Logger::info << "JointSearch was called as follow:" << endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << endl << endl;
}

void GecoArguments::printSummary() {
  string boolStr[2] = {string("OFF"), string("ON")};
  Logger::info << "Parameters summary: " << endl;
  Logger::info << "Gene tree: " << geneTree << endl;
  Logger::info << "Alignment: " << alignment << endl; 
  Logger::info << "Species tree: " << speciesTree << endl;
  Logger::info << "Gene species map: " << geneSpeciesMap << endl;
  Logger::info << "Strategy: " << strategy << endl;
  Logger::info << "Reconciliation model: " << reconciliationModel << endl;
  Logger::info << "Reconciliation opt: " << reconciliationOpt << endl;
  Logger::info << "Libpll model: " << libpllModel << endl;
  Logger::info << "Prefix: " << output << endl;
  Logger::info << "Check mode: " << boolStr[check] << endl;
  Logger::info << "Unrooted gene tree: " << boolStr[!rootedGeneTree] << endl;
  Logger::info << "MPI Ranks: " << ParallelContext::getSize() << endl;
  Logger::info << endl;
}
