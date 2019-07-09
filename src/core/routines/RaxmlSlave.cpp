#include "RaxmlSlave.hpp"
#include <string>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/ParallelOfstream.hpp>
#include <parallelization/ParallelContext.hpp>
#include <likelihoods/LibpllEvaluation.hpp>

int RaxmlSlave::runRaxmlOptimization(int argc, char** argv, void* comm)
{
  assert(argc == 8);
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

