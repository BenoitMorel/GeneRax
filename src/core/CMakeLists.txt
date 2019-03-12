
set(jointsearch-core_SOURCES ParallelContext.cpp
  Scenario.cpp
  PerCoreGeneTrees.cpp
  IO/Logger.cpp
  IO/GeneSpeciesMapping.cpp
  IO/FamiliesFileParser.cpp
  IO/ParallelOfstream.cpp
  likelihoods/LibpllEvaluation.cpp
  likelihoods/ReconciliationEvaluation.cpp
  likelihoods/reconciliation_models/AbstractReconciliationModel.cpp
  likelihoods/reconciliation_models/DatedDLModel.cpp
  likelihoods/reconciliation_models/UndatedDLModel.cpp
  likelihoods/reconciliation_models/UndatedDTLModel.cpp
  optimizers/DTLOptimizer.cpp
  treeSearch/JointTree.cpp
  treeSearch/Moves.cpp
  treeSearch/Rollbacks.cpp
  treeSearch/SearchUtils.cpp
  treeSearch/SPRSearch.cpp
  )

add_library(jointsearch-core ${jointsearch-core_SOURCES})

target_include_directories(jointsearch-core
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )


