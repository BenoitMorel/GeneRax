
set(jointsearch-core_SOURCES 
  IO/Families.cpp
  IO/Logger.cpp
  IO/GeneSpeciesMapping.cpp
  IO/FamiliesFileParser.cpp
  IO/LibpllParsers.cpp
  IO/Model.cpp
  IO/ParallelOfstream.cpp
  IO/ReconciliationWriter.cpp
  likelihoods/LibpllEvaluation.cpp
  likelihoods/ReconciliationEvaluation.cpp
  maths/Random.cpp
  NJ/MiniNJ.cpp
  NJ/Cherry.cpp
  NJ/CherryPro.cpp
  optimizers/DTLOptimizer.cpp
  optimizers/PerFamilyDTLOptimizer.cpp
  optimizers/SpeciesTreeOptimizer.cpp
  parallelization/ParallelContext.cpp
  parallelization/PerCoreGeneTrees.cpp
  parallelization/Scheduler.cpp
  routines/scheduled_routines/GeneRaxSlave.cpp
  routines/scheduled_routines/GeneRaxMaster.cpp
  routines/scheduled_routines/RaxmlMaster.cpp
  routines/scheduled_routines/RaxmlSlave.cpp
  routines/Routines.cpp
  routines/SlavesMain.cpp
  search/Moves.cpp
  search/Rollbacks.cpp
  search/SearchUtils.cpp
  search/SPRSearch.cpp
  trees/Clade.cpp
  trees/PLLUnrootedTree.cpp
  trees/PLLRootedTree.cpp
  trees/PLLTreeInfo.cpp
  trees/JointTree.cpp
  trees/SpeciesTree.cpp
  trees/TreeDuplicatesFinder.cpp
  util/Scenario.cpp
  )

add_library(jointsearch-core STATIC ${jointsearch-core_SOURCES})

target_include_directories(jointsearch-core
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )


