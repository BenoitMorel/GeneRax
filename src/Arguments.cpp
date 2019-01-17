#include "Arguments.hpp"
#include <Logger.hpp>
#include <ParallelContext.hpp>

int Arguments::argc = 0;
char ** Arguments::argv = 0;
string Arguments::geneTree;
string Arguments::alignment;
string Arguments::speciesTree;
string Arguments::geneSpeciesMap;
string Arguments::strategy("SPR");
Arguments::ReconciliationModel Arguments::reconciliationModel(UndatedDL);
string Arguments::output("jointSearch");
bool Arguments::check = false;
bool Arguments::rootedGeneTree = false;
bool Arguments::noFelsensteinLikelihood = false;  

Arguments::ReconciliationModel getModelFromString(const string &modelStr)
{
  if (modelStr == "UndatedDL") {
    return Arguments::UndatedDL;
  } else if (modelStr == "DatedDL") {
    return Arguments::DatedDL;
  } else if (modelStr == "UndatedDTL") {
    return Arguments::UndatedDTL;
  } else {
    return Arguments::InvalidModel;
  }
}

string getStringFromModel(Arguments::ReconciliationModel model) {
  if (model == Arguments::UndatedDL) {
    return "UndatedDL";
  } else if (model == Arguments::UndatedDTL) {
    return "UndatedDTL";
  } else if (model == Arguments::DatedDL) {
    return "DatedDL";
  } else {
    return "InvalidModel";
  }
}

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
    } else if (arg == "-m" || arg == "--map") {
      geneSpeciesMap = string(argv[++i]);
    } else if (arg == "--strategy") {
      strategy = string(argv[++i]);
    } else if (arg == "--reconciliation-model") {
      reconciliationModel = getModelFromString(string(argv[++i]));
    } else if (arg == "-p" || arg == "--prefix") {
      output = string(argv[++i]);
    } else if (arg == "--check") {
      check = true;
    } else if (arg == "--rooted-gene-tree") {
      rootedGeneTree = true;
    } else if (arg == "--no-felsenstein-likelihood") {
      noFelsensteinLikelihood = true;
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
  if (!geneSpeciesMap.size()) {
    Logger::error << "You need to provide a gene species map file." << endl;
    ok = false;
  }
  if (reconciliationModel == Arguments::InvalidModel) {
    Logger::error << "Invalid reconciliation model." << endl;
    ok = false;
  }
  if (!ok) {
    Logger::error << "Aborting." << endl;
    exit(1);
  }
  assertFileExists(geneTree);
  assertFileExists(speciesTree);
  assertFileExists(geneSpeciesMap);
  assertFileExists(alignment);
}

void Arguments::printHelp() {
  Logger::info << "-h, --help" << endl;
  Logger::info << "-g, --gene-tree <GENE TREE>" << endl;
  Logger::info << "-a, --alignment <ALIGNMENT>" << endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << endl;
  Logger::info << "--strategy <STRATEGY>  {EVAL, SPR}" << endl;
  Logger::info << "--reconciliation-model <reconciliationModel>  {UndatedDL, UndatedDTL, DatedDL}" << endl;
  Logger::info << "-t, --threads <THREADS NUMBER>" << endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << endl;
  Logger::info << "--check" << endl;
  Logger::info << "--rooted-gene-tree" << endl;
  Logger::info << "--no-felsenstein-likelihood" << endl;
  Logger::info << endl;

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
  Logger::info << "Gene species map: " << geneSpeciesMap << endl;
  Logger::info << "Strategy: " << strategy << endl;
  Logger::info << "Reconciliation model: " << getStringFromModel(reconciliationModel) << endl;
  Logger::info << "Prefix: " << output << endl;
  Logger::info << "Check mode: " << boolStr[check] << endl;
  Logger::info << "Rooted gene tree: " << boolStr[rootedGeneTree] << endl;
  Logger::info << "Ignoring felsenstein likelihood: " << boolStr[noFelsensteinLikelihood] << endl;
  Logger::info << "MPI Ranks: " << ParallelContext::getSize() << endl;
  Logger::info << endl;
}
