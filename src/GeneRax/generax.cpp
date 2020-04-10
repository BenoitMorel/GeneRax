#include "GeneRaxArguments.hpp"
#include "GeneRaxCore.hpp"
#include "GeneRaxInstance.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <routines/SlavesMain.hpp>




int generax_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm); 
  Logger::init();
  Logger::timed << "GeneRax v1.1.0" << std::endl; 
  GeneRaxInstance instance(argc, argv);
  GeneRaxCore::initInstance(instance);
  GeneRaxCore::initRandomGeneTrees(instance);
  GeneRaxCore::printStats(instance);
  GeneRaxCore::speciesTreeSearch(instance);
  GeneRaxCore::rerootSpeciesTree(instance);
  GeneRaxCore::geneTreeJointSearch(instance);
  GeneRaxCore::reconcile(instance);
  GeneRaxCore::terminate(instance);
  
  Logger::close();
  ParallelContext::finalize();
  return 0;
}


/**
 * GeneRax can call itself for executing routines from the scheduler. This function
 * decides whether we are in the main called by the user (generax_main)
 * or if we are called by the scheduler to execute some 
 * intermediate routine (static_scheduled_main)
*/
int internal_main(int argc, char** argv, void* comm)
{
  if (SlavesMain::isSlave(argc, argv)) {
    int slaveComm = -1; 
    return static_scheduled_main(argc, argv, &slaveComm);
  } else {
    return generax_main(argc, argv, comm);
  }
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

