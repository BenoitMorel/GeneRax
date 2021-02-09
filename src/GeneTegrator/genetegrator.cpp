#include "GeneTegratorArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "GTSpeciesTreeOptimizer.hpp"



void run(const GeneTegratorArguments &args)
{
  Logger::info << "Mkdir " << args.output << std::endl;
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  GTSpeciesTreeOptimizer speciesTreeOptimizer(
      args.speciesTree,
      families,
      args.output);
  auto ll = speciesTreeOptimizer.computeRecLikelihood();
  Logger::info << "total ll=" << ll << std::endl;
  ll = speciesTreeOptimizer.sprSearch(1);
  
  Logger::info << "total ll=" << ll << std::endl;
}

int internal_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm); 
  Logger::init();
  Logger::timed << "GeneTegrator v0.0.0" << std::endl; 
  GeneTegratorArguments args(argc, argv); 
  run(args);
  Logger::close();
  ParallelContext::finalize();
  return 0;
}
int main(int argc, char** argv)
{
#ifdef WITH_MPI
  return internal_main(argc, argv, 0);
#else
  int noMPIComm = -1;
  return internal_main(argc, argv, &noMPIComm);
#endif
}


