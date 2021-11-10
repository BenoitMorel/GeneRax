
#include "GeneRaxCore.hpp"
#include "GeneRaxInstance.hpp"
#include <parallelization/ParallelContext.hpp>
#include <branchlengths/ReconciliationBLEstimator.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <random>
#include <limits>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <maths/Parameters.hpp>
#include <IO/FileSystem.hpp>
#include <IO/ParallelOfstream.hpp>
#include <maths/Random.hpp>
#include <NJ/MiniNJ.hpp>
#include <NJ/Cherry.hpp>
#include <NJ/CherryPro.hpp>
#include <parallelization/Scheduler.hpp>
#include <routines/Routines.hpp>
#include <optimizers/SpeciesTreeOptimizer.hpp>
#include <trees/SpeciesTree.hpp>
#include <support/ICCalculator.hpp>
#include <util/Paths.hpp>

static void initStartingSpeciesTree(GeneRaxInstance &instance)
{
  instance.speciesTree = Paths::getSpeciesTreeFile(
      instance.args.output, 
      "starting_species_tree.newick");
  std::unique_ptr<PLLRootedTree> speciesTree(nullptr);
  if (instance.args.speciesTreeAlgorithm == SpeciesTreeAlgorithm::User) {
    unsigned int canRead = 1;
    if (ParallelContext::getRank() == 0) {
      try {
        SpeciesTree reader(instance.args.speciesTree);
      } catch (const std::exception &e) {
        Logger::info << "Error while trying to parse the species tree:" << std::endl;
        Logger::info << e.what() << std::endl;
        canRead = 0;
      }
    }
    ParallelContext::broadcastUInt(0, canRead);
    if (!canRead) {
      ParallelContext::abort(153);
    }
    // add labels to internal nodes
    LibpllParsers::labelRootedTree(instance.args.speciesTree, instance.speciesTree);
  } else {
    Routines::computeInitialSpeciesTree(instance.currentFamilies,
        instance.args.output,
        instance.args.speciesTreeAlgorithm)->save(instance.speciesTree);

  }
  ParallelContext::barrier();
  if (ParallelContext::getRank() == 0) {
    SpeciesTree copy(instance.speciesTree); 
    instance.speciesTree = Paths::getSpeciesTreeFile(instance.args.output, "inferred_species_tree.newick");
    copy.getTree().save(instance.speciesTree);
  }
  ParallelContext::barrier();
}

void GeneRaxCore::initInstance(GeneRaxInstance &instance) 
{
  Random::setSeed(static_cast<unsigned int>(instance.args.seed));
  FileSystem::mkdir(instance.args.output, true);
  Logger::initFileOutput(FileSystem::joinPaths(instance.args.output, "generax"));
  // assert twice, before of a bug I had at the 
  // second rand() call with openmpi
  assert(ParallelContext::isRandConsistent());
  assert(ParallelContext::isRandConsistent());
  instance.args.printCommand();
  instance.args.printSummary();
  instance.initialFamilies = FamiliesFileParser::parseFamiliesFile(instance.args.families);
  initFolders(instance);
  bool needAlignments = instance.args.strategy != GeneSearchStrategy::SKIP
    && instance.args.strategy != GeneSearchStrategy::RECONCILE;
  if (instance.args.filterFamilies) {
    Logger::timed << "Filtering invalid families..." << std::endl;
    bool checkSpeciesTree = (instance.args.speciesTreeAlgorithm 
        == SpeciesTreeAlgorithm::User);
    Family::filterFamilies(instance.initialFamilies, instance.args.speciesTree, needAlignments, checkSpeciesTree);
    if (!instance.initialFamilies.size()) {
      Logger::info << "[Error] No valid families! Aborting GeneRax" << std::endl;
      ParallelContext::abort(10);
    }
  }
  instance.currentFamilies = instance.initialFamilies;
  initFolders(instance);
  instance.modelParameters = ModelParameters(instance.rates, 
      instance.currentFamilies.size(),
      instance.getRecModelInfo());
  Logger::info << "End of instance initialization" << std::endl;
}

void GeneRaxCore::initRandomGeneTrees(GeneRaxInstance &instance)
{
  assert(ParallelContext::isRandConsistent());
  instance.currentFamilies = instance.initialFamilies;
  bool randoms = Routines::createRandomTrees(instance.args.output, instance.currentFamilies); 
  if (randoms) {
    initialGeneTreeSearch(instance);
  }
}

void GeneRaxCore::initSpeciesTree(GeneRaxInstance &instance) 
{
  Logger::timed << "Starting species tree initialization..." << std::endl;
  initStartingSpeciesTree(instance);
  Logger::timed << "End of species tree initialization" << std::endl;
  if (instance.args.filterFamilies) {
    Logger::timed << "Filtering invalid families based on the starting species tree..." << std::endl;
    bool needAlignments = instance.args.strategy != GeneSearchStrategy::SKIP
      && instance.args.strategy != GeneSearchStrategy::RECONCILE;
    Family::filterFamilies(instance.currentFamilies, instance.speciesTree, needAlignments, true);
  }
  if (!instance.currentFamilies.size()) {
    Logger::info << "[Error] No valid families! Aborting GeneRax" << std::endl;
    ParallelContext::abort(10);
  }
  
}
  
void GeneRaxCore::generateFakeAlignments(GeneRaxInstance &instance)
{
  if (!instance.args.generateFakeAlignments) {
    return;
  }
  Logger::timed << "Generating fake alignments" << std::endl;
  std::string fakeDir = FileSystem::joinPaths(instance.args.output, "fake_msas");
  FileSystem::mkdir(fakeDir, true);
  ParallelContext::barrier();
  PerCoreGeneTrees perCoreTrees(instance.currentFamilies);
  std::unordered_set<std::string> coreFamilies;
  for (const auto &perCoreTree: perCoreTrees.getTrees()) {
    coreFamilies.insert(perCoreTree.name);
  }
  for (auto &family: instance.currentFamilies) {
    family.alignmentFile = FileSystem::joinPaths(fakeDir, family.name);
    family.libpllModel = "GTR";
    if (coreFamilies.find(family.name) != coreFamilies.end()) {
      // generate the MSA
      PLLUnrootedTree tree(family.startingGeneTree);
      std::ofstream os(family.alignmentFile);
      for (auto leaf: tree.getLeaves()) {
        os << ">" << leaf->label << std::endl << "ACGT" << std::endl;
      }
    }
  }
}
  
void GeneRaxCore::printStats(GeneRaxInstance &instance)
{
  if (instance.args.filterFamilies) {
    std::string coverageFile = FileSystem::joinPaths(instance.args.output,
      std::string("perSpeciesCoverage.txt"));
    std::string fractionMissingFile = FileSystem::joinPaths(instance.args.output,
      std::string("fractionMissing.txt"));
    Logger::timed << "Gathering statistics about the families..." << std::endl;
    Family::printStats(instance.currentFamilies, 
      instance.speciesTree,
      coverageFile,
      fractionMissingFile);
  }
}


static void speciesTreeSearchAux(GeneRaxInstance &instance, int samples)
{
  Families saveFamilies = instance.currentFamilies;
  if (samples > 0) {
    auto rng = std::default_random_engine {};
    std::shuffle(instance.currentFamilies.begin(), instance.currentFamilies.end(), rng);
    instance.currentFamilies.resize(samples);
  }

  ParallelContext::barrier();
  Parameters startingRates = instance.getUserParameters();
  SpeciesTreeSearchParams searchParams;
  searchParams.sprRadius = instance.args.speciesSPRRadius;
  searchParams.rootSmallRadius = instance.args.speciesSmallRootRadius;
  searchParams.rootBigRadius = instance.args.speciesBigRootRadius;
  SpeciesTreeOptimizer speciesTreeOptimizer(instance.speciesTree, 
      instance.currentFamilies, 
      instance.getRecModelInfo(), 
      startingRates, 
      instance.args.userDTLRates, 
      instance.args.output, 
      searchParams);
  if (instance.args.speciesSPRRadius > 0) {
    Logger::info << std::endl;
    Logger::timed << "Start optimizing the species tree with fixed gene trees (on " 
      << instance.currentFamilies.size() << " families " << std::endl;
  }
  speciesTreeOptimizer.optimize(instance.args.speciesStrategy);
  instance.totalRecLL = speciesTreeOptimizer.getReconciliationLikelihood();
  instance.speciesTree = speciesTreeOptimizer.saveCurrentSpeciesTreeId();
  instance.totalRecLL = speciesTreeOptimizer.getReconciliationLikelihood();
  Logger::info << std::endl;
  Logger::timed << "[Species search] End of optimizing the species tree" << std::endl;
  instance.speciesTree = speciesTreeOptimizer.saveCurrentSpeciesTreeId();

  instance.currentFamilies = saveFamilies;
  ParallelContext::barrier();
}

void GeneRaxCore::speciesTreeSearch(GeneRaxInstance &instance)
{
  assert(ParallelContext::isRandConsistent());
  if (instance.args.speciesStrategy == SpeciesSearchStrategy::SKIP) {
    return;
  }
  Logger::info << "Saving tree to " << instance.speciesTree << std::endl;
  speciesTreeSearchAux(instance, -1);
}


void GeneRaxCore::geneTreeJointSearch(GeneRaxInstance &instance)
{
  assert(ParallelContext::isRandConsistent());
  if (instance.args.strategy == GeneSearchStrategy::SKIP ||
      instance.args.strategy == GeneSearchStrategy::RECONCILE) {
    return;
  }
  for (unsigned int i = 1; i <= instance.args.recRadius; ++i) { 
    bool enableLibpll = false;
    bool perSpeciesDTLRates = false;
    optimizeRatesAndGeneTrees(instance, perSpeciesDTLRates, enableLibpll, i);
  }
  for (unsigned int i = 1; i <= instance.args.maxSPRRadius; ++i) {
    bool enableLibpll = true;
    bool perSpeciesDTLRates = instance.args.perSpeciesDTLRates && (i >= instance.args.maxSPRRadius - 1); // only apply per-species optimization at the two last rounds
    optimizeRatesAndGeneTrees(instance, perSpeciesDTLRates, enableLibpll, i);
  }
  ModelParameters modelRates(instance.rates,
      1,
      instance.getRecModelInfo());
  instance.readModelParameters(modelRates);
  instance.modelParameters = modelRates;
}


void GeneRaxCore::reconcile(GeneRaxInstance &instance)
{
  assert(ParallelContext::isRandConsistent());
  if (instance.args.reconcile || instance.args.reconciliationSamples > 0) {
    if (instance.args.strategy == GeneSearchStrategy::RECONCILE) {
      Logger::timed << "Optimizing DTL rates before the reconciliation..." << std::endl;
      // we haven't optimized the DTL rates yet, so we do it now
      if (!instance.args.perFamilyDTLRates) {
        Routines::optimizeRates(instance.args.userDTLRates, 
          instance.speciesTree, 
          instance.recModelInfo,
          instance.currentFamilies, 
          instance.args.perSpeciesDTLRates, 
          instance.rates, 
          instance.elapsedRates);
        instance.totalRecLL = instance.rates.getScore();
      } else {
        long elapsed = 0;
        bool enableLibpll = false;
        unsigned int sprRadius = 0;
        Routines::optimizeGeneTrees(
            instance.currentFamilies, 
            instance.recModelInfo,
            instance.rates, 
            instance.args.output, 
            "results", 
            instance.args.execPath, 
            instance.speciesTree, 
            RecOpt::Grid,
            false,
            instance.args.supportThreshold, 
            instance.args.recWeight, 
            true, 
            enableLibpll, 
            sprRadius, 
            instance.currentIteration++, 
            ParallelContext::allowSchedulerSplitImplementation(), 
            elapsed);
        double temp = 0.0;
        Routines::gatherLikelihoods(instance.currentFamilies, temp, instance.totalRecLL);
      }
    }
        
    Logger::timed << "Reconciling gene trees with the species tree..." << std::endl;
    bool optimizeRates = false;
    Routines::inferReconciliation(instance.speciesTree, 
        instance.currentFamilies, 
        instance.modelParameters, 
        instance.args.output, 
        instance.args.reconcile,
        instance.args.reconciliationSamples, 
        optimizeRates);
    if (instance.args.buildSuperMatrix) {
      std::string outputSuperMatrixAll = FileSystem::joinPaths(
          instance.args.output, "superMatrixAll.fasta");
      Routines::computeSuperMatrixFromOrthoGroups(instance.speciesTree,
        instance.currentFamilies,
        instance.args.output, 
        outputSuperMatrixAll,
        false,
        true);
    }
  }
}
  
void GeneRaxCore::terminate(GeneRaxInstance &instance)
{
  assert(ParallelContext::isRandConsistent());
  Logger::timed << "Terminating the instance.." << std::endl;
  ParallelOfstream os(FileSystem::joinPaths(instance.args.output, "stats.txt"));
  os << "JointLL: " << instance.totalLibpllLL + instance.totalRecLL << std::endl;
  os << "LibpllLL: " << instance.totalLibpllLL << std::endl;
  os << "RecLL: " << instance.totalRecLL;
  Logger::info << std::endl;
  auto &rates = instance.rates;
  if (!instance.args.perFamilyDTLRates) {
    if (rates.dimensions() == 2) {
      Logger::timed << "DT rates: D=" << rates[0] << " L= " << rates[1] << std::endl;
    } else if (instance.rates.dimensions() == 3) {
      Logger::timed<< "DTL rates: D=" << rates[0] << " L= " << rates[1] << " T=" << rates[2] << std::endl;
    }
  }
  Logger::timed << "Reconciliation likelihood: " << instance.totalRecLL << std::endl;
  if (instance.totalLibpllLL) {
    Logger::timed << "Phylogenetic likelihood: " << instance.totalLibpllLL << std::endl;
    Logger::timed << "Joint likelihood: " << instance.totalLibpllLL + instance.totalRecLL << std::endl;
  }
#ifdef PRINT_TIMES
  if (instance.elapsedRaxml) {
    Logger::timed << "Initial time spent on optimizing random trees: " << instance.elapsedRaxml << "s" << std::endl;
  }
  Logger::timed << "Time spent on optimizing rates: " << instance.elapsedRates << "s" << std::endl;
  Logger::timed << "Time spent on optimizing gene trees: " << instance.elapsedSPR << "s" << std::endl;
#endif
  Logger::timed << "Results directory: " << instance.args.output << std::endl;
  Logger::timed << "End of GeneRax execution" << std::endl;
}


void GeneRaxCore::initFolders(GeneRaxInstance &instance) 
{
  assert(ParallelContext::isRandConsistent());
  std::string results = FileSystem::joinPaths(instance.args.output, "results");
  FileSystem::mkdir(results, true);
  for (auto &family: instance.currentFamilies) {
    FileSystem::mkdir(FileSystem::joinPaths(results, family.name), true);
  }
  for (auto dir: Paths::getDirectoriesToCreate(instance.args.output)) {
    FileSystem::mkdir(dir, true);
  }
}

void GeneRaxCore::initialGeneTreeSearch(GeneRaxInstance &instance)
{
  assert(ParallelContext::isRandConsistent());
  Logger::info << std::endl;
  Logger::timed << "[Initialization] Initial optimization of the starting random gene trees" << std::endl;
  Logger::timed << "[Initialization] All the families will first be optimized with sequences only" << std::endl;
  Logger::mute();
  Routines::runRaxmlOptimization(instance.currentFamilies, instance.args.output, 
      instance.args.execPath, instance.currentIteration++, 
      ParallelContext::allowSchedulerSplitImplementation(), instance.elapsedRaxml);
  Logger::unmute();
  Routines::gatherLikelihoods(instance.currentFamilies, instance.totalLibpllLL, instance.totalRecLL);
  Logger::timed << "[Initialization] Finished optimizing some of the gene trees" << std::endl;
  Logger::info << std::endl;
}

void GeneRaxCore::optimizeRatesAndGeneTrees(GeneRaxInstance &instance,
    bool perSpeciesDTLRates,
    bool enableLibpll,
    unsigned int sprRadius)
{
  assert(ParallelContext::isRandConsistent());
  long elapsed = 0;
  if (!instance.args.perFamilyDTLRates) {
    Logger::timed << "Reconciliation rates optimization... " << std::endl;
    Routines::optimizeRates(instance.args.userDTLRates, 
        instance.speciesTree, 
        instance.recModelInfo,
        instance.currentFamilies, 
        perSpeciesDTLRates, 
        instance.rates, 
        instance.elapsedRates);
    if (perSpeciesDTLRates) {
      auto rateFile = FileSystem::joinPaths(instance.args.output, "per_species_rates.txt");
      Routines::exportPerSpeciesRates(instance.speciesTree,
          instance.rates,
          instance.recModelInfo,
          rateFile);
    }
    if (!instance.args.perFamilyDTLRates && !instance.args.perSpeciesDTLRates) {
      auto paramNames = Enums::parameterNames(instance.recModelInfo.model);
      Logger::info << "\t";
      for (unsigned int i = 0; i < paramNames.size(); ++i) {
        Logger::info << paramNames[i] << "=" << instance.rates[i] << ", ";
      }
      Logger::info << "RecLL= " << instance.rates.getScore() << std::endl;
    } else {
      Logger::info << "\tRecLL=" << instance.rates.getScore() << std::endl;
    }
    Logger::info << std::endl;
  }
  std::string additionalMsg;
  if (instance.args.perFamilyDTLRates) {
    additionalMsg = std::string("reconciliation rates and ");
  }
  Logger::timed << "Optimizing " + additionalMsg + "gene trees with radius=" << sprRadius << "... " << std::endl; 
  const bool enableRecLL = true;
  Routines::optimizeGeneTrees(instance.currentFamilies, 
      instance.recModelInfo, 
      instance.rates, 
      instance.args.output, 
      "results", 
      instance.args.execPath, 
      instance.speciesTree, 
      RecOpt::Grid, 
      instance.args.madRooting,
      instance.args.supportThreshold, 
      instance.args.recWeight, 
      enableRecLL, 
      enableLibpll, 
      sprRadius, 
      instance.currentIteration++, 
      ParallelContext::allowSchedulerSplitImplementation(), 
      elapsed);
  instance.elapsedSPR += elapsed;
  Routines::gatherLikelihoods(instance.currentFamilies, instance.totalLibpllLL, instance.totalRecLL);
  Logger::info << "\tJointLL=" << instance.totalLibpllLL + instance.totalRecLL 
    << " RecLL=" << instance.totalRecLL << " LibpllLL=" << instance.totalLibpllLL << std::endl;
  Logger::info << std::endl;
}
  
void GeneRaxCore::speciesTreeBLEstimation(GeneRaxInstance &instance)
{
  if (!instance.args.estimateSpeciesBranchLenghts) {
    return;
  }
  ParallelContext::barrier();
  ReconciliationBLEstimator::estimate(
      instance.speciesTree,
      instance.currentFamilies,
      instance.modelParameters);
  ParallelContext::barrier();
}

void GeneRaxCore::speciesTreeSupportEstimation(GeneRaxInstance &instance)
{
 
  ParallelContext::barrier();
  if (instance.args.quartetSupport) { 
    ICCalculator calculator(instance.speciesTree,
        instance.currentFamilies,
        instance.args.eqpicRadius,
        true);
    auto qpicOutput = Paths::getSpeciesTreeFile(instance.args.output,
        "species_tree_qpic.newick");
    auto eqpicOutput = Paths::getSpeciesTreeFile(instance.args.output,
        "species_tree_eqpic.newick");
    auto supportOutput = Paths::getSpeciesTreeFile(instance.args.output,
        "species_tree_quartet_support.newick");
    auto supportTripletOutput = Paths::getSpeciesTreeFile(
        instance.args.output,
        "species_tree_quartet_support_triplet.newick");
    calculator.exportScores(qpicOutput, 
        eqpicOutput, 
        supportOutput,
        supportTripletOutput);
    ParallelContext::barrier();
    FileSystem::copy(eqpicOutput, instance.speciesTree, true);
     
  }
  ParallelContext::barrier();


  if (instance.args.quartetSupportAllQuartets) { 
    Logger::timed 
      << "Start estimating species tree support values" 
      << std::endl;
    ICCalculator calculator(instance.speciesTree,
        instance.currentFamilies,
        instance.args.eqpicRadius,
        false);
    auto qpicOutput = Paths::getSpeciesTreeFile(instance.args.output,
        "species_tree_qpic_allquartets.newick");
    auto eqpicOutput = Paths::getSpeciesTreeFile(instance.args.output,
        "species_tree_eqpic_allquartets.newick");
    auto supportOutput = Paths::getSpeciesTreeFile(instance.args.output,
        "species_tree_quartet_support_allquartets.newick");
    auto supportTripletOutput = Paths::getSpeciesTreeFile(
        instance.args.output,
        "species_tree_quartet_support_triplet_allquartets.newick");
    calculator.exportScores(qpicOutput, 
        eqpicOutput, 
        supportOutput,
        supportTripletOutput);
    Logger::timed 
      << "Finished estimating species tree support values" 
      << std::endl;
  }
  ParallelContext::barrier();

}



