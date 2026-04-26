#include "GeneRaxArguments.hpp"
#include "GeneRaxCore.hpp"
#include "GeneRaxInstance.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <routines/SlavesMain.hpp>
#include <util/Checkpoint.hpp>
#include <IO/FileSystem.hpp>
#include <util/Paths.hpp>


int generax_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm); 
  Logger::init();
  Logger::timed << "GeneRax 2.1.3" << std::endl; 
  std::vector<double> plop;
  plop.assign(1, 1.0);
  GeneRaxInstance instance(argc, argv);
  GeneRaxCore::initInstance(instance);
  Checkpoint checkpoint(instance.args.outputPath, instance.args.checkpoint);
  bool resuming = checkpoint.isEnabled() && checkpoint.has("phase");
  int completedPhase = resuming ? checkpoint.getInt("phase") : -1;
  if (resuming) {
    Logger::timed << "[Checkpoint] Resuming from phase " << completedPhase << std::endl;
  }
  if (completedPhase < 1) {
    GeneRaxCore::initRandomGeneTrees(instance);
    GeneRaxCore::initSpeciesTree(instance);
    if (checkpoint.isEnabled()) {
      // Save the starting species tree so we can fall back to it on
      // resume when no iteration-level checkpoint exists yet
      auto startingTree = Paths::getSpeciesTreeFile(
          instance.args.outputPath, "checkpoint_starting_species_tree.newick");
      FileSystem::copy(instance.speciesTree, startingTree, true);
      checkpoint.save("phase", 1);
    }
  } else {
    // Resuming: pick the right species tree depending on whether
    // an iteration-level checkpoint exists
    auto iterationTree = Paths::getSpeciesTreeFile(
        instance.args.outputPath, "checkpoint_species_tree.newick");
    auto startingTree = Paths::getSpeciesTreeFile(
        instance.args.outputPath, "checkpoint_starting_species_tree.newick");
    if (checkpoint.has("hybrid_iteration") && FileSystem::exists(iterationTree)) {
      // Resume from the tree saved at the last completed search iteration
      instance.speciesTree = iterationTree;
      Logger::timed << "[Checkpoint] Loaded species tree from iteration checkpoint"
        << std::endl;
    } else if (FileSystem::exists(startingTree)) {
      // No iteration checkpoint yet — restart search from the starting tree
      instance.speciesTree = startingTree;
      Logger::timed << "[Checkpoint] Loaded starting species tree from checkpoint"
        << std::endl;
    } else {
      Logger::timed << "[Checkpoint] Warning: no checkpointed species tree found, reinitializing" << std::endl;
      GeneRaxCore::initRandomGeneTrees(instance);
      GeneRaxCore::initSpeciesTree(instance);
    }
  }
  GeneRaxCore::generateFakeAlignments(instance);
  GeneRaxCore::printStats(instance);
  if (completedPhase < 2) {
    GeneRaxCore::speciesTreeSearch(instance);
    if (checkpoint.isEnabled()) {
      checkpoint.save("phase", 2);
    }
  } else {
    Logger::timed << "[Checkpoint] Skipping species tree search (already completed)" << std::endl;
  }
  GeneRaxCore::geneTreeJointSearch(instance);
  GeneRaxCore::reconcile(instance);
  GeneRaxCore::speciesTreeBLEstimation(instance);
  GeneRaxCore::speciesTreeSupportEstimation(instance);
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

