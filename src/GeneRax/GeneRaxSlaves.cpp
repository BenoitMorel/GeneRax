#include "GeneRaxArguments.hpp"
#include <ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <limits>
#include <PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/DTLRates.hpp>
#include <treeSearch/JointTree.hpp>
#include <treeSearch/SPRSearch.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <sstream>

void getTreeStrings(const string &filename, vector<string> &treeStrings) 
{
  string geneTreeString;
  if (filename == "__random__" || filename.size() == 0) {
    treeStrings.push_back("__random__");
    return;
  }
  ifstream treeStream(filename);
  while(getline(treeStream, geneTreeString)) {
    geneTreeString.erase(remove(geneTreeString.begin(), geneTreeString.end(), '\n'), geneTreeString.end());
    if (geneTreeString.empty()) {
      continue;
    }
    treeStrings.push_back(geneTreeString);
  }
}


void optimizeGeneTreesSlave(const string &startingGeneTreeFile,
    const string &mappingFile,
    const string &alignmentFile,
    const string &speciesTreeFile,
    const string libpllModel, 
    RecModel recModel,
    RecOpt recOpt,
    bool rootedGeneTree,
    double dupRate,
    double lossRate, 
    double transferRate,
    bool enableRec,
    int sprRadius,
    const string &outputGeneTree,
    const string &outputStats) 
{
  Logger::timed << "Starting optimizing gene tree" << endl;
  vector<string> geneTreeStrings;
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
      false, //check
      false,
      dupRate,
      lossRate,
      transferRate
      );
  jointTree->enableReconciliation(enableRec);
  Logger::info << "Taxa number: " << jointTree->getGeneTaxaNumber() << endl;
  jointTree->optimizeParameters(true, false); // only optimize felsenstein likelihood
  double bestLoglk = jointTree->computeJointLoglk();
  jointTree->printLoglk();
  Logger::info << "Initial ll = " << bestLoglk << endl;
  if (sprRadius > 0) {
    while(SPRSearch::applySPRRound(*jointTree, sprRadius, bestLoglk, true)) {} 
  }
  Logger::info << "Final ll = " << bestLoglk << endl;
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
  Logger::timed << "End of optimizing gene tree" << endl;
  ParallelContext::barrier();
}


int optimizeGeneTreesMain(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  if (argc != 17) {
    Logger::error << "Invalid number of parameters in generax_optimize_gene_trees: " << argc << endl;
    return 1;
  }
  Logger::timed << "Starting optimizeGeneTreesSlave" << endl;
  int i = 2;
  string startingGeneTreeFile(argv[i++]);
  string mappingFile(argv[i++]);
  string alignmentFile(argv[i++]);
  string speciesTreeFile(argv[i++]);
  string libpllModel(argv[i++]);
  Logger::info << "LibpllModel " << libpllModel << endl;
  RecModel recModel = RecModel(atoi(argv[i++])); 
  RecOpt recOpt = RecOpt(atoi(argv[i++])); 
  bool rootedGeneTree = bool(atoi(argv[i++]));
  double dupRate = double(atof(argv[i++]));
  double lossRate = double(atof(argv[i++]));
  double transferRate = double(atof(argv[i++]));
  bool enableRec = bool(atoi(argv[i++]));
  int sprRadius = atoi(argv[i++]);
  string outputGeneTree(argv[i++]);
  string outputStats(argv[i++]);
  optimizeGeneTreesSlave(startingGeneTreeFile,
      mappingFile,
      alignmentFile,
      speciesTreeFile,
      libpllModel,
      recModel,
      recOpt,
      rootedGeneTree,
      dupRate,
      lossRate,
      transferRate,
      enableRec,
      sprRadius,
      outputGeneTree,
      outputStats);
  ParallelContext::finalize();
  Logger::timed << "End of optimizeGeneTreesSlave" << endl;
  return 0;
}

int raxmlLightMain(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  int i = 2;
  string startingGeneTreeFile(argv[i++]);
  string alignmentFile(argv[i++]);
  string libpllModel(argv[i++]);
  string outputGeneTree(argv[i++]);
  string outputLibpllModel(argv[i++]);
  string outputStats(argv[i++]);
  LibpllAlignmentInfo info;
  info.alignmentFilename = alignmentFile;
  info.model = libpllModel;
  Logger::info << startingGeneTreeFile << endl;
  auto evaluation = LibpllEvaluation::buildFromFile(startingGeneTreeFile,
      info);
  double ll = 0.0;
  evaluation->optimizeAllParameters(10.0);
  evaluation->raxmlSPRRounds(1, 5, 0);
  ll = evaluation->raxmlSPRRounds(1, 15, 0);
  Logger::info << evaluation->computeLikelihood(false) << endl; 
  ll = evaluation->raxmlSPRRounds(1, 15, 0);
  Logger::info << evaluation->computeLikelihood(false) << endl;
  ll = evaluation->raxmlSPRRounds(1, 20, 0);
  Logger::info << evaluation->computeLikelihood(false) << endl;
  ll = evaluation->raxmlSPRRounds(1, 25, 0);
  Logger::info << evaluation->computeLikelihood(false) << endl;
  Logger::info << "optimize all parameters" << endl;
  evaluation->optimizeAllParameters(3.0);
  
  
  ll = evaluation->raxmlSPRRounds(1, 5, 1);
  Logger::info << evaluation->computeLikelihood(false) << endl;
  ll = evaluation->raxmlSPRRounds(1, 10, 1);
  Logger::info << evaluation->computeLikelihood(false) << endl;
  ll = evaluation->raxmlSPRRounds(1, 15, 1);
  
  evaluation->optimizeAllParameters(1.0);
  
  
  Logger::info << evaluation->computeLikelihood(false) << endl;
  /*
  ll = evaluation->raxmlSPRRounds(1, 15, 1);
  ll = evaluation->raxmlSPRRounds(1, 20, 1);
  */
  ll = evaluation->raxmlSPRRounds(1, 25, 1);
  Logger::info << evaluation->computeLikelihood(false) << endl;
  
  ParallelOfstream stats(outputStats);
  stats << ll <<  " 0.0";
  stats.close();
  string modelStr = evaluation->getModelStr();
  
  ParallelOfstream modelWriter(outputLibpllModel);
  modelWriter << modelStr << endl;
  modelWriter.close();
  LibpllParsers::saveUtree(evaluation->getTreeInfo()->root, outputGeneTree);  
  ParallelContext::finalize();
  return 0;
}


extern "C" int static_scheduled_main(int argc, char** argv, void* comm)
{
  Logger::enableLogFile(false);
  string key(argv[1]);
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





