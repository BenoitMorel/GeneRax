
set(jointsearch-core_SOURCES 
  IO/Logger.cpp
  IO/GeneSpeciesMapping.cpp
  IO/FamiliesFileParser.cpp
  IO/LibpllParsers.cpp
  IO/Model.cpp
  IO/ParallelOfstream.cpp
  likelihoods/LibpllEvaluation.cpp
  likelihoods/ReconciliationEvaluation.cpp
  likelihoods/SubtreeRepeatsCache.cpp
  optimizers/DTLOptimizer.cpp
  optimizers/PerFamilyDTLOptimizer.cpp
  optimizers/SpeciesTreeOptimizer.cpp
  parallelization/ParallelContext.cpp
  parallelization/Scheduler.cpp
  routines/GeneRaxSlaves.cpp
  routines/RaxmlMaster.cpp
  routines/RaxmlSlave.cpp
  routines/Routines.cpp
  routines/GeneTreeSearchMaster.cpp
  search/Moves.cpp
  search/Rollbacks.cpp
  search/SearchUtils.cpp
  search/SPRSearch.cpp
  trees/JointTree.cpp
  trees/PerCoreGeneTrees.cpp
  trees/SpeciesTree.cpp
  util/Scenario.cpp
  )

add_library(jointsearch-core STATIC ${jointsearch-core_SOURCES})

target_include_directories(jointsearch-core
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )


