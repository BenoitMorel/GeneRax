/**
 * This is a simpler reimplementation of the raxml-ng search algorithm with autodetect, fast and slow SPR roudns
 */


#include "RaxmlSlave.hpp"
#include <string>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/ParallelOfstream.hpp>
#include <parallelization/ParallelContext.hpp>
#include <likelihoods/LibpllEvaluation.hpp>

static void optimizeParameters(LibpllEvaluation &evaluation, double radius) 
{
  double initialLL = evaluation.computeLikelihood(false);
  Logger::timed << "[" << initialLL << "]" << " Optimize parametres (" << radius << ")" << std::endl;
  evaluation.optimizeAllParameters(10.0);
}

static void optimizeBranches(LibpllEvaluation &evaluation, double radius) 
{
  evaluation.optimizeBranches(radius);
}

static bool optimizeTopology(LibpllEvaluation &evaluation, 
    unsigned int radiusMin, 
    unsigned int radiusMax, 
    unsigned int thorough, 
    unsigned int toKeep, 
    double cutoff) 
{
  
  double initialLL = evaluation.computeLikelihood(false);
  Logger::timed << "["  << initialLL << "] " 
    << (thorough ? "SLOW" : "FAST") << " SPR Round (radius=" << radiusMax << ")" << std::endl;
  double ll = evaluation.raxmlSPRRounds(radiusMin, radiusMax, thorough, toKeep, cutoff);
  optimizeBranches(evaluation, 1.0);
  return ll > initialLL + 0.1;
}

int RaxmlSlave::runRaxmlOptimization(int argc, char** argv, void* comm)
{
  assert(argc == 8);
  ParallelContext::init(comm);
  Logger::init();
  int i = 2;
  std::string startingGeneTreeFile(argv[i++]);
  std::string alignmentFile(argv[i++]);
  std::string libpllModel(argv[i++]);
  std::string outputGeneTree(argv[i++]);
  std::string outputLibpllModel(argv[i++]);
  std::string outputStats(argv[i++]);
  Logger::info << startingGeneTreeFile << std::endl;
  LibpllEvaluation evaluation(startingGeneTreeFile, true, alignmentFile, libpllModel);
  Logger::timed << "LL = " << evaluation.computeLikelihood(false) << std::endl; 
  optimizeBranches(evaluation, 1.0);
  optimizeParameters(evaluation, 10.0);
  unsigned int radiusMin = 1;
  unsigned int radiusMax = 5;
  unsigned int radiusLimit = std::min(22u, 
      static_cast<unsigned int>(evaluation.getTreeInfo()->tip_count - 3) );
  unsigned int toKeep = 0;
  unsigned int thorough = 0;
  double cutoff = 0.0;
  while (radiusMax < radiusLimit) {
    if (optimizeTopology(evaluation, radiusMin, radiusMax, thorough, toKeep, cutoff)) {
      radiusMin += 5;
      radiusMax += 5;
    } else {
      radiusMin -= 5;
      radiusMax -= 5;
      break;
    }
  }
  optimizeParameters(evaluation, 3.0);
  radiusMin = 1;
  toKeep = 20;
  cutoff = 1.0;
  while (optimizeTopology(evaluation, radiusMin, radiusMax, 0, toKeep, cutoff)) {
  }
  thorough = 1;
  optimizeParameters(evaluation, 1.0);
  for (radiusMax = 5; radiusMax <= radiusLimit; radiusMax += 5) {
    if (optimizeTopology(evaluation, 1, radiusMax, thorough, toKeep, cutoff)) {
      radiusMax = 5;
    }
  }
  ParallelOfstream stats(outputStats);
  stats << evaluation.computeLikelihood(false) <<  " 0.0";
  stats.close();
  std::string modelStr = evaluation.getModelStr();
  ParallelOfstream modelWriter(outputLibpllModel);
  modelWriter << modelStr << std::endl;
  modelWriter.close();
  LibpllParsers::saveUtree(evaluation.getTreeInfo()->root, outputGeneTree);  
  ParallelContext::finalize();
  return 0;
}

