#include "GeneRaxSlave.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <limits>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/Parameters.hpp>
#include <trees/JointTree.hpp>
#include <search/SPRSearch.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <sstream>
#include <routines/scheduled_routines/RaxmlSlave.hpp>

static void getTreeStrings(const std::string &filename, std::vector<std::string> &treeStrings) 
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


static void optimizeGeneTreesSlave(const std::string &startingGeneTreeFile,
    const std::string &mappingFile,
    const std::string &alignmentFile,
    const std::string &speciesTreeFile,
    const std::string &libpllModel, 
    const std::string &ratesFile,
    const RecModelInfo &recModelInfo,
    RecOpt recOpt,
    bool madRooting,
    double supportThreshold,
    double recWeight,
    bool enableRec,
    bool enableLibpll,
    int sprRadius,
    const std::string &outputGeneTree,
    const std::string &outputStats) 
{
  auto start = std::chrono::high_resolution_clock::now();
  Logger::timed << "Starting optimizing gene tree" << std::endl;
  Logger::info << "Number of ranks " << ParallelContext::getSize() << std::endl;
  std::vector<std::string> geneTreeStrings;
  getTreeStrings(startingGeneTreeFile, geneTreeStrings);
  assert(geneTreeStrings.size() == 1);
  Parameters ratesVector(ratesFile);
  auto jointTree = std::make_unique<JointTree>(geneTreeStrings[0],
      alignmentFile,
      speciesTreeFile,
      mappingFile,
      libpllModel,
      recModelInfo,
      recOpt,
      madRooting,
      supportThreshold,
      recWeight,
      false, //check
      recModelInfo.perFamilyRates,
      ratesVector
      );
  jointTree->enableReconciliation(enableRec);
  jointTree->enableLibpll(enableLibpll);
  Logger::info << "Taxa number: " << jointTree->getGeneTaxaNumber() << std::endl;
  jointTree->optimizeParameters(true,  enableRec);
  double bestLoglk = jointTree->computeJointLoglk();
  jointTree->printLoglk();
  Logger::info << "Initial ll = " << bestLoglk << std::endl;
  if (sprRadius > 0) {
    while(SPRSearch::applySPRRound(*jointTree, sprRadius, true)) {} 
  }
  jointTree->printLoglk();
  if (outputGeneTree.size() && ParallelContext::getRank() == 0) {
    Logger::info << "Saving tree in " <<outputGeneTree << std::endl;
    jointTree->save(outputGeneTree, false);
    Logger::info << "Finished saving  " <<outputGeneTree << std::endl;

  }
  if (outputStats.size()) {
    ParallelOfstream stats(outputStats);
    double libpllLL = jointTree->computeLibpllLoglk ();
    double recLL = jointTree->computeReconciliationLoglk();
    stats << libpllLL << " " << recLL << std::endl;
    stats << "Reconciliation rates = ";
    for (auto rate: jointTree->getRatesVector().getVector()) {
      stats << rate << " ";
    }
    stats.close();
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  auto seconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
  Logger::timed << "End of optimizing gene tree after " << seconds << "s" << std::endl;
  ParallelContext::barrier();
}

static std::string getArg(const std::string &str)
{
  return (str == "NONE" ? std::string() : str);
}

int GeneRaxSlave::optimizeGeneTreesMain(int argc, char** argv, void* comm)
{
  assert(argc == 17 + RecModelInfo::getArgc());
  ParallelContext::init(comm);
  Logger::timed << "Starting optimizeGeneTreesSlave" << std::endl;
  int i = 2;
  std::string startingGeneTreeFile(argv[i++]);
  std::string mappingFile(getArg(argv[i++]));
  std::string alignmentFile(argv[i++]);
  if (alignmentFile == "NOALIGNMENT") {
    alignmentFile = "";
  }
  std::string speciesTreeFile(argv[i++]);
  std::string libpllModel(argv[i++]);
  std::string ratesFile(argv[i++]);
  Logger::info << "LibpllModel " << libpllModel << std::endl;
  RecModelInfo recModelInfo;
  recModelInfo.readFromArgv(argv, i);
  RecOpt recOpt = RecOpt(atoi(argv[i++])); 
  
  //recModelInfo.perFamilyRates = bool(atoi(argv[i++]));
  //recModel rootedGeneTree = bool(atoi(argv[i++]));
  double supportThreshold = double(atof(argv[i++]));
  double recWeight = double(atof(argv[i++]));
  bool enableRec = bool(atoi(argv[i++]));
  bool enableLibpll = bool(atoi(argv[i++]));
  int sprRadius = atoi(argv[i++]);
  std::string outputGeneTree(argv[i++]);
  std::string outputStats(argv[i++]);
  bool madRooting = bool(atoi(argv[i++]));
  optimizeGeneTreesSlave(startingGeneTreeFile,
      mappingFile,
      alignmentFile,
      speciesTreeFile,
      libpllModel,
      ratesFile,
      recModelInfo,
      recOpt,
      madRooting,
      supportThreshold,
      recWeight,
      enableRec,
      enableLibpll,
      sprRadius,
      outputGeneTree,
      outputStats);
  ParallelContext::finalize();
  Logger::timed << "End of optimizeGeneTreesSlave" << std::endl;
  return 0;
}


