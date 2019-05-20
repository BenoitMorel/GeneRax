
set(jointsearch-core_SOURCES ParallelContext.cpp
  Scenario.cpp
  IO/Logger.cpp
  IO/GeneSpeciesMapping.cpp
  IO/FamiliesFileParser.cpp
  IO/LibpllParsers.cpp
  IO/Model.cpp
  IO/ParallelOfstream.cpp
  likelihoods/LibpllEvaluation.cpp
  likelihoods/ReconciliationEvaluation.cpp
  likelihoods/reconciliation_models/AbstractReconciliationModel.cpp
  likelihoods/reconciliation_models/UndatedDLModel.cpp
  likelihoods/reconciliation_models/UndatedDTLModel.cpp
  optimizers/DTLOptimizer.cpp
  trees/JointTree.cpp
  trees/PerCoreGeneTrees.cpp
  search/Moves.cpp
  search/Rollbacks.cpp
  search/SearchUtils.cpp
  search/SPRSearch.cpp
  )

add_library(jointsearch-core STATIC ${jointsearch-core_SOURCES})

target_include_directories(jointsearch-core
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )


