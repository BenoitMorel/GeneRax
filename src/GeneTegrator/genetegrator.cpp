#include "GeneTegratorArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "UndatedDLMultiModel.hpp"
#include <trees/PLLRootedTree.hpp>

void run(const GeneTegratorArguments &args)
{
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  PLLRootedTree speciesTree(args.speciesTree);
  double sumLL = 0.0;
  for (auto &family: families) {
    ConditionalClades ccp(family.startingGeneTree);    
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    UndatedDLMultiModel<double> model(speciesTree,
        mapping,
        ccp);
    double ll = model.computeLogLikelihood();
    Logger::info << "ll=" << ll;
    sumLL += ll;
  }
  Logger::info << "total ll = " << sumLL << std::endl;
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


