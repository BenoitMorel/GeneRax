
set(jointsearch-core_SOURCES 
  branchlengths/ReconciliationBLEstimator.cpp
  ccp/ConditionalClades.cpp
  ccp/RootedSpeciesSplitScore.cpp
  ccp/SpeciesSplits.cpp
  ccp/UnrootedSpeciesSplitScore.cpp
  IO/NewickParserCommon.cpp
  IO/RootedNewickParser.cpp
  IO/Families.cpp
  IO/HighwayCandidateParser.cpp
  IO/IO.cpp
  IO/Logger.cpp
  IO/GeneSpeciesMapping.cpp
  IO/FamiliesFileParser.cpp
  IO/LibpllParsers.cpp
  IO/Model.cpp
  IO/ParallelOfstream.cpp
  IO/ReconciliationWriter.cpp
  likelihoods/LibpllEvaluation.cpp
  likelihoods/ReconciliationEvaluation.cpp
  likelihoods/reconciliation_models/BaseReconciliationModel.cpp
  maths/bitvector.cpp
  maths/Random.cpp
  DistanceMethods/Astrid.cpp
  DistanceMethods/Asteroid.cpp
  DistanceMethods/MiniNJ.cpp
  DistanceMethods/Cherry.cpp
  DistanceMethods/CherryPro.cpp
  DistanceMethods/NeighborJoining.cpp
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
  search/SpeciesSearchCommon.cpp
  search/SpeciesRootSearch.cpp
  search/SpeciesSPRSearch.cpp
  search/SpeciesTransferSearch.cpp
  search/SPRSearch.cpp
  search/UNNISearch.cpp
  search/DatedSpeciesTreeSearch.cpp
  support/ICCalculator.cpp
  trees/Clade.cpp
  trees/DSTagger.cpp
  trees/DatedTree.cpp
  trees/PLLUnrootedTree.cpp
  trees/PLLRootedTree.cpp
  trees/PLLTreeInfo.cpp
  trees/PolytomySolver.cpp
  trees/PolyTree.cpp
  trees/JointTree.cpp
  trees/SpeciesTree.cpp
  util/Scenario.cpp
  )

add_library(jointsearch-core STATIC ${jointsearch-core_SOURCES})

target_include_directories(jointsearch-core
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  )


