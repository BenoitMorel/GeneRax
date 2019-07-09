#include "GeneRaxSlaves.hpp"
#include "GeneRaxArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <limits>
#include <trees/PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/DTLRates.hpp>
#include <trees/JointTree.hpp>
#include <search/SPRSearch.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <sstream>
#include <routines/RaxmlSlave.hpp>

void getTreeStrings(const std::string &filename, std::vector<std::string> &treeStrings) 
{
  std::string geneTreeString;
  if (filename == "__random__" || filename.size() == 0) {
    treeStrings.push_back("__random__");
    return;
  }
  std::ifstream treeStream(filename);
  while(getline(treeStream, geneTreeString)) {
    geneTreeString.erase(remove(geneTreeString.begin(), geneTreeString.end(), '\n'), geneTreeString.end());
    if (geneTreeString.empty()) {
      continue;
    }
    treeStrings.push_back(geneTreeString);
  }
}


void optimizeGeneTreesSlave(const std::string &startingGeneTreeFile,
    const std::string &mappingFile,
    const std::string &alignmentFile,
    const std::string &speciesTreeFile,
    const std::string &libpllModel, 
    const std::string &ratesFile,
    RecModel recModel,
    RecOpt recOpt,
    bool rootedGeneTree,
    double recWeight,
    bool enableRec,
    int sprRadius,
    const std::string &outputGeneTree,
    const std::string &outputStats) 
{
  Logger::timed << "Starting optimizing gene tree" << std::endl;
  std::vector<std::string> geneTreeStrings;
  getTreeStrings(startingGeneTreeFile, geneTreeStrings);
  assert(geneTreeStrings.size() == 1);
  DTLRatesVector ratesVector(ratesFile);
  auto jointTree = std::make_shared<JointTree>(geneTreeStrings[0],
      alignmentFile,
      speciesTreeFile,
      mappingFile,
      libpllModel,
      recModel,
      recOpt,
      rootedGeneTree,
      recWeight,
      false, //check
      false, // optimize DTL
      ratesVector
      );
  jointTree->enableReconciliation(enableRec);
  Logger::info << "Taxa number: " << jointTree->getGeneTaxaNumber() << std::endl;
  jointTree->optimizeParameters(true, false); // only optimize felsenstein likelihood
  double bestLoglk = jointTree->computeJointLoglk();
  jointTree->printLoglk();
  Logger::info << "Initial ll = " << bestLoglk << std::endl;
  if (sprRadius > 0) {
    while(SPRSearch::applySPRRound(*jointTree, sprRadius, bestLoglk, true)) {} 
  }
  Logger::info << "Final ll = " << bestLoglk << std::endl;
  if (outputGeneTree.size()) {
    jointTree->save(outputGeneTree, false);
  }
  if (outputStats.size()) {
    ParallelOfstream stats(outputStats);
    double libpllLL = jointTree->computeLibpllLoglk ();
    double recLL = jointTree->computeReconciliationLoglk();
    stats << libpllLL << " " << recLL;
    stats.close();
  }
  Logger::timed << "End of optimizing gene tree" << std::endl;
  ParallelContext::barrier();
}

std::string getArg(const std::string &str)
{
  return (str == "NONE" ? std::string() : str);
}

int optimizeGeneTreesMain(int argc, char** argv, void* comm)
{
  assert(argc == 16);
  ParallelContext::init(comm);
  Logger::timed << "Starting optimizeGeneTreesSlave" << std::endl;
  int i = 2;
  std::string startingGeneTreeFile(argv[i++]);
  std::string mappingFile(getArg(argv[i++]));
  std::string alignmentFile(argv[i++]);
  std::string speciesTreeFile(argv[i++]);
  std::string libpllModel(argv[i++]);
  std::string ratesFile(argv[i++]);
  Logger::info << "LibpllModel " << libpllModel << std::endl;
  RecModel recModel = RecModel(atoi(argv[i++])); 
  RecOpt recOpt = RecOpt(atoi(argv[i++])); 
  bool rootedGeneTree = bool(atoi(argv[i++]));
  double recWeight = double(atof(argv[i++]));
  bool enableRec = bool(atoi(argv[i++]));
  int sprRadius = atoi(argv[i++]);
  std::string outputGeneTree(argv[i++]);
  std::string outputStats(argv[i++]);
  optimizeGeneTreesSlave(startingGeneTreeFile,
      mappingFile,
      alignmentFile,
      speciesTreeFile,
      libpllModel,
      ratesFile,
      recModel,
      recOpt,
      rootedGeneTree,
      recWeight,
      enableRec,
      sprRadius,
      outputGeneTree,
      outputStats);
  ParallelContext::finalize();
  Logger::timed << "End of optimizeGeneTreesSlave" << std::endl;
  return 0;
}


bool GeneRaxSlaves::is_slave(int argc, char** argv)
{
  if (argc < 2) {
    return false;
  }
  std::string key(argv[1]);
  return key == "optimizeGeneTrees" || key == "raxmlLight";
}

extern "C" int static_scheduled_main(int argc, char** argv, void* comm)
{
  Logger::enableLogFile(false);
  std::string key(argv[1]);
  int res = 1;
  if (key == "optimizeGeneTrees") {
    res = optimizeGeneTreesMain(argc, argv, comm);   
  } else if (key == "raxmlLight") {
    res = RaxmlSlave::runRaxmlOptimization(argc, argv, comm);
  } else {
    assert(0);
  }
  Logger::enableLogFile(true);
  return res;
}

