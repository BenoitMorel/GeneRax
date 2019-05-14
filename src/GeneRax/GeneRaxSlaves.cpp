#include "GeneRaxArguments.hpp"
#include <ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <limits>
#include <trees/PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/DTLRates.hpp>
#include <trees/JointTree.hpp>
#include <treeSearch/SPRSearch.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <sstream>

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
    const std::string libpllModel, 
    RecModel recModel,
    RecOpt recOpt,
    bool rootedGeneTree,
    double recWeight,
    double dupRate,
    double lossRate, 
    double transferRate,
    bool enableRec,
    int sprRadius,
    const std::string &outputGeneTree,
    const std::string &outputStats) 
{
  Logger::timed << "Starting optimizing gene tree" << std::endl;
  std::vector<std::string> geneTreeStrings;
  getTreeStrings(startingGeneTreeFile, geneTreeStrings);
  assert(geneTreeStrings.size() == 1);
  auto jointTree = make_shared<JointTree>(geneTreeStrings[0],
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
      dupRate,
      lossRate,
      transferRate
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


int optimizeGeneTreesMain(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  if (argc != 18) {
    Logger::error << "Invalid number of parameters in generax_optimize_gene_trees: " << argc << std::endl;
    return 1;
  }
  Logger::timed << "Starting optimizeGeneTreesSlave" << std::endl;
  int i = 2;
  std::string startingGeneTreeFile(argv[i++]);
  std::string mappingFile(argv[i++]);
  std::string alignmentFile(argv[i++]);
  std::string speciesTreeFile(argv[i++]);
  std::string libpllModel(argv[i++]);
  Logger::info << "LibpllModel " << libpllModel << std::endl;
  RecModel recModel = RecModel(atoi(argv[i++])); 
  RecOpt recOpt = RecOpt(atoi(argv[i++])); 
  bool rootedGeneTree = bool(atoi(argv[i++]));
  double recWeight = double(atof(argv[i++]));
  double dupRate = double(atof(argv[i++]));
  double lossRate = double(atof(argv[i++]));
  double transferRate = double(atof(argv[i++]));
  bool enableRec = bool(atoi(argv[i++]));
  int sprRadius = atoi(argv[i++]);
  std::string outputGeneTree(argv[i++]);
  std::string outputStats(argv[i++]);
  optimizeGeneTreesSlave(startingGeneTreeFile,
      mappingFile,
      alignmentFile,
      speciesTreeFile,
      libpllModel,
      recModel,
      recOpt,
      rootedGeneTree,
      recWeight,
      dupRate,
      lossRate,
      transferRate,
      enableRec,
      sprRadius,
      outputGeneTree,
      outputStats);
  ParallelContext::finalize();
  Logger::timed << "End of optimizeGeneTreesSlave" << std::endl;
  return 0;
}

int raxmlLightMain(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  int i = 2;
  std::string startingGeneTreeFile(argv[i++]);
  std::string alignmentFile(argv[i++]);
  std::string libpllModel(argv[i++]);
  std::string outputGeneTree(argv[i++]);
  std::string outputLibpllModel(argv[i++]);
  std::string outputStats(argv[i++]);
  LibpllAlignmentInfo info;
  info.alignmentFilename = alignmentFile;
  info.model = libpllModel;
  Logger::info << startingGeneTreeFile << std::endl;
  auto evaluation = LibpllEvaluation::buildFromFile(startingGeneTreeFile,
      info);
  double ll = 0.0;
  evaluation->optimizeAllParameters(10.0);
  evaluation->raxmlSPRRounds(1, 5, 0);
  ll = evaluation->raxmlSPRRounds(1, 15, 0);
  Logger::info << evaluation->computeLikelihood(false) << std::endl; 
  ll = evaluation->raxmlSPRRounds(1, 15, 0);
  Logger::info << evaluation->computeLikelihood(false) << std::endl;
  ll = evaluation->raxmlSPRRounds(1, 20, 0);
  Logger::info << evaluation->computeLikelihood(false) << std::endl;
  ll = evaluation->raxmlSPRRounds(1, 25, 0);
  Logger::info << evaluation->computeLikelihood(false) << std::endl;
  Logger::info << "optimize all parameters" << std::endl;
  evaluation->optimizeAllParameters(3.0);
  
  
  ll = evaluation->raxmlSPRRounds(1, 5, 1);
  Logger::info << evaluation->computeLikelihood(false) << std::endl;
  ll = evaluation->raxmlSPRRounds(1, 10, 1);
  Logger::info << evaluation->computeLikelihood(false) << std::endl;
  ll = evaluation->raxmlSPRRounds(1, 15, 1);
  
  evaluation->optimizeAllParameters(1.0);
  
  
  Logger::info << evaluation->computeLikelihood(false) << std::endl;
  /*
  ll = evaluation->raxmlSPRRounds(1, 15, 1);
  ll = evaluation->raxmlSPRRounds(1, 20, 1);
  */
  ll = evaluation->raxmlSPRRounds(1, 25, 1);
  Logger::info << evaluation->computeLikelihood(false) << std::endl;
  
  ParallelOfstream stats(outputStats);
  stats << ll <<  " 0.0";
  stats.close();
  std::string modelStr = evaluation->getModelStr();
  
  ParallelOfstream modelWriter(outputLibpllModel);
  modelWriter << modelStr << std::endl;
  modelWriter.close();
  LibpllParsers::saveUtree(evaluation->getTreeInfo()->root, outputGeneTree);  
  ParallelContext::finalize();
  return 0;
}


extern "C" int static_scheduled_main(int argc, char** argv, void* comm)
{
  Logger::enableLogFile(false);
  std::string key(argv[1]);
  int res = 1;
  if (key == "optimizeGeneTrees") {
    res = optimizeGeneTreesMain(argc, argv, comm);   
  } else if (key == "raxmlLight") {
    res = raxmlLightMain(argc, argv, comm);
  } else {
    assert(0);
  }
  Logger::enableLogFile(true);
  return res;
}


#ifdef SLAVE_EXEC
int main(int argc, char** argv) {
  int comm = -1;
  return static_scheduled_main(argc, argv, &comm);
}
#endif
