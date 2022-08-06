/**
 * This is a simpler reimplementation of the raxml-ng search algorithm with autodetect, fast and slow SPR roudns
 */


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
    unsigned int radius, 
    unsigned int thorough, 
    unsigned int toKeep, 
    double cutoff) 
{
  
  double initialLL = evaluation.computeLikelihood(false);
  Logger::timed << "["  << initialLL << "] " 
    << (thorough ? "SLOW" : "FAST") << " SPR Round (radius=" << radius << ")" << std::endl;
  double ll = evaluation.raxmlSPRRounds(1, radius, thorough, toKeep, cutoff);
  Logger::info << "ll=" << ll << std::endl;
  return ll > initialLL + 0.1;
}

int lightSearch(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  Logger::init();
  if (argc != 7) {
    Logger::info << "Syntax error" << std::endl;
    Logger::info << "lightsearch starting_gene_tree alignment model radius maxiterations outputtree" 
      << std::endl;
    return 1;
  }
  int i = 1;
  std::string startingGeneTreeFile(argv[i++]);
  std::string alignmentFile(argv[i++]);
  std::string libpllModel(argv[i++]);
  int radius = atoi(argv[i++]);
  int maxIterations = atoi(argv[i++]);
  std::string outputTree(argv[i++]);
  
  
  
  LibpllEvaluation evaluation(startingGeneTreeFile, true, alignmentFile, libpllModel);
  Logger::timed << "LL = " << evaluation.computeLikelihood(false) << std::endl; 
  optimizeBranches(evaluation, 1.0);
  optimizeParameters(evaluation, 10.0);
  unsigned int toKeep = 0;
  unsigned int thorough = 0;
  double cutoff = 0.0;
  for (int it = 0; it < maxIterations; ++it) {
    if (!optimizeTopology(evaluation, radius, thorough, toKeep, cutoff)) {
      break;
    }
  }
  LibpllParsers::saveUtree(evaluation.getTreeInfo()->root, outputTree);  
  ParallelContext::finalize();
  return 0;
}

int main(int argc, char** argv) {
  return lightSearch(argc, argv, nullptr);
}
