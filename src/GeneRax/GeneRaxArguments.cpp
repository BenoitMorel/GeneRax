
#include "GeneRaxArguments.hpp"
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>
#include <algorithm>
#include <vector>
#include <util/Constants.hpp>

GeneRaxArguments::GeneRaxArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  geneSearchStrategy(GeneSearchStrategy::SPR),
  reconciliationModelStr("UndatedDTL"),
  outputPath("GeneRax"),
  perFamilyDTLRates(false),
  rootedGeneTree(true),
  madRooting(false),
  pruneSpeciesTree(false),
  supportThreshold(-1.0),
  recRadius(0),
  perSpeciesDTLRates(false),
  userDTLRates(false),
  transferConstraint(TransferConstaint::PARENTS),
  noDup(false),
  dupRate(0.2),
  lossRate(0.2),
  transferRate(0.2),
  gammaCategories(1),
  originationStrategy(OriginationStrategy::UNIFORM),
  reconcile(true),
  buildSuperMatrix(false),
  reconciliationSampleNumber(0),
  maxSPRRadius(5),
  recWeight(1.0), 
  seed(123),
  filterFamilies(true),
  exec(iargv[0]),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  speciesStrategy(SpeciesSearchStrategy::SKIP),
  constrainSpeciesSearch(false),
  rerootSpeciesTree(false),
  estimateSpeciesBranchLenghts(false), 
  speciesSPRRadius(DEFAULT_SPECIES_SPR_RADIUS),
  speciesSmallRootRadius(DEFAULT_SPECIES_SMALL_ROOT_RADIUS),
  speciesBigRootRadius(DEFAULT_SPECIES_BIG_ROOT_RADIUS),
  minGeneBranchLength(0.000001),
  quartetSupport(false),
  quartetSupportAllQuartets(false),
  eqpicRadius(std::numeric_limits<int>::max()),
  generateFakeAlignments(false)
{
  if (argc == 1) {
    printHelp();
    ParallelContext::abort(0);
  }
  init();
}

void GeneRaxArguments::init() {
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    /**
     * General parameters
     */
    if (arg == "-h" || arg == "--help") {
      printHelp();
      ParallelContext::abort(0);
    } else if (arg == "-f" || arg == "--families") {
      familyFilePath = std::string(argv[++i]);
    } else if (arg == "-s" || arg == "--species-tree") {
      speciesTree = std::string(argv[++i]);
      speciesTreeAlgorithm = Enums::strToSpeciesTree(speciesTree);
    } else if (arg == "--strategy") {
      geneSearchStrategy = ArgumentsHelper::strToStrategy(std::string(argv[++i]));
      if (geneSearchStrategy == GeneSearchStrategy::EVAL) {
        recRadius = maxSPRRadius = 0;
      }
    } else if (arg == "-p" || arg == "--prefix") {
      outputPath = std::string(argv[++i]);
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--skip-family-filtering") {
      filterFamilies = false;
    /**
     *  Model parameters
     */
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModelStr = std::string(argv[++i]);
    } else if (arg == "--per-family-rates") {
      perFamilyDTLRates = true;
    } else if (arg == "--unrooted-gene-tree") {
      rootedGeneTree = false;
    } else if (arg == "--mad-rooting") {
      madRooting = true;
    } else if (arg == "--prune-species-tree") {
      pruneSpeciesTree = true;
    } else if (arg == "--rec-radius") {
      recRadius = static_cast<unsigned int>(atoi(argv[++i]));
    } else if (arg == "--per-species-rates") {
      perSpeciesDTLRates = true;
    } else if (arg == "--dtl-rates-opt") {
      if (ArgumentsHelper::strToRecOpt(argv[++i]) == 
          RecOpt::None) {
        userDTLRates = true;
      }
    } else if (arg == "--transfer-constraint") {
      transferConstraint = ArgumentsHelper::strToTransferConstraint(std::string(argv[++i]));
    } else if (arg == "--no-dup") {
      dupRate = 0.0;
      noDup = true;
    } else if (arg == "--dup-rate") {
      dupRate = atof(argv[++i]);
    } else if (arg == "--loss-rate") {
      lossRate = atof(argv[++i]);
    } else if (arg == "--transfer-rate") {
      transferRate = atof(argv[++i]);
    } else if (arg == "--rec-weight") {
      recWeight = atof(argv[++i]);
    } else if (arg == "--reconciliation-samples") {
      reconciliationSampleNumber = static_cast<unsigned int>(atoi(argv[++i]));
    /**
     * Gene tree correction
     */
    } else if (arg == "--geneSearchStrategy") {
      geneSearchStrategy = ArgumentsHelper::strToStrategy(std::string(argv[++i]));
      if (geneSearchStrategy == GeneSearchStrategy::EVAL) {
        recRadius = maxSPRRadius = 0;
      }
    } else if (arg == "--reconcile") {
      reconcile = true;
    } else if (arg == "--do-not-reconcile") {
      reconcile = false;
    } else if (arg == "--support-threshold") {
      supportThreshold = static_cast<double>(atof(argv[++i]));
    } else if (arg == "--max-spr-radius") {
      maxSPRRadius = static_cast<unsigned int>(atoi(argv[++i]));
    /**
     *  Species tree inference
     */
    } else if (arg == "--si-strategy") {
      speciesStrategy = ArgumentsHelper::strToSpeciesSearchStrategy(std::string(argv[++i]));
    } else if (arg == "--si-spr-radius") {
      speciesSPRRadius = static_cast<unsigned int>(atoi(argv[++i]));
    } else if (arg == "--si-small-root-radius") {
      speciesSmallRootRadius = static_cast<unsigned int>(atoi(argv[++i]));
    } else if (arg == "--si-big-root-radius") {
      speciesBigRootRadius = static_cast<unsigned int>(atoi(argv[++i]));
    } else if (arg == "--si-constrained-search") {
      constrainSpeciesSearch = true;
    } else if (arg == "--si-estimate-bl") {
      estimateSpeciesBranchLenghts = true;
    } else if (arg == "--si-quartet-support") {
      quartetSupport = true; 
    } else if (arg == "--si-eqpic-radius") {
      eqpicRadius = atoi(argv[++i]); 
    /**
     *  Experimental and debug
     */
    } else if (arg == "--fraction-missing") {
      fractionMissingFile = argv[++i];
    } else if (arg == "--build-supermatrix") {
      reconcile = true;
      buildSuperMatrix = true;
    } else if (arg == "--gene-min-bl") {
      minGeneBranchLength = atof(argv[++i]);
    } else if (arg == "--quartet-support-all-quartets") {
      quartetSupportAllQuartets = true; 
    } else if (arg == "--reconciliation-samples") {
      reconciliationSampleNumber = static_cast<unsigned int>(atoi(argv[++i]));
    } else if (arg == "--generate-fake-alignments") {
      generateFakeAlignments = true;
    } else {
      Logger::info << "Unrecognized argument " << arg << std::endl;
      Logger::info << "Aborting" << std::endl;
      ParallelContext::abort(1);
    }
  }
  execPath = std::string(argv[0]);
  checkInputs();
}

static void assertFileExists(const std::string &file) 
{
  std::ifstream f(file);
  if (!f) {
    Logger::info << "File " << file << " does not exist. Aborting." << std::endl;
    ParallelContext::abort(0);
  }
}


void GeneRaxArguments::checkInputs() {
  bool ok = true;
  assertFileExists(familyFilePath); 
  if (!speciesTree.size() && speciesStrategy != SpeciesSearchStrategy::SKIP) {
    Logger::info << "[Error] You need to provide a species tree or to optimize it." << std::endl;
    ok = false;
  }
  if (userDTLRates && perSpeciesDTLRates) {
    Logger::info << "[Error] You cannot specify the rates when using per-species DTL rates" << std::endl;
    ok = false;
  }
  if (userDTLRates && (dupRate < 0.0 || lossRate < 0.0)) {
    Logger::info << "[Error] You specified at least one of the duplication and loss rates, but not both of them." << std::endl;
    ok = false;
  }
  if (perSpeciesDTLRates && perFamilyDTLRates) {
    Logger::info << "[Error] You cannot use per-family and per-species rates at the same time" << std::endl;
    ok = false;
  }
  if (!ArgumentsHelper::isValidRecModel(reconciliationModelStr)) {
    Logger::info << "[Error] Invalid reconciliation model string " << reconciliationModelStr << std::endl;
    ok = false;
  }
  if (!ok) {
    Logger::info << "Aborting." << std::endl;
    ParallelContext::abort(1);
  }
  if (speciesTreeAlgorithm == SpeciesTreeAlgorithm::User) {
    assertFileExists(speciesTree);
  }
}

void GeneRaxArguments::printHelp() {
  Logger::info << "-h, --help" << std::endl;
  Logger::info << "-f, --families <FAMILIES_INFORMATION>" << std::endl;
  Logger::info << "-s, --species-tree <SPECIES TREE>" << std::endl;
  Logger::info << "--geneSearchStrategy <STRATEGY>  {EVAL, SPR}" << std::endl;
  Logger::info << "-r --rec-model <reconciliationModel>  {UndatedDL, UndatedDTL, Auto}" << std::endl;
  Logger::info << "-p, --prefix <OUTPUT PREFIX>" << std::endl;
  Logger::info << "--unrooted-gene-tree" << std::endl;
  Logger::info << "--support-threshold <threshold>" << std::endl;
  Logger::info << "--per-family-rates" << std::endl;
  Logger::info << "--per-species-rates" << std::endl;
  Logger::info << "--dup-rate <duplication rate>" << std::endl;
  Logger::info << "--loss-rate <loss rate>" << std::endl;
  Logger::info << "--transfer-rate <transfer rate>" << std::endl;
  Logger::info << "--max-spr-radius <max SPR radius>" << std::endl;
  Logger::info << "--rec-weight <reconciliation likelihood weight>" << std::endl;
  Logger::info << "--do-not-reconcile" << std::endl;
  Logger::info << "--reconciliation-samples <number of samples>" << std::endl;
  Logger::info << "--seed <seed>" << std::endl;
  Logger::info << "Please find more information on the GeneRax github wiki" << std::endl;
  Logger::info << std::endl;

}

void GeneRaxArguments::printCommand() {
  Logger::info << "GeneRax was called as follow:" << std::endl;
  for (int i = 0; i < argc; ++i) {
    Logger::info << argv[i] << " ";
  }
  Logger::info << std::endl << std::endl;
}

void GeneRaxArguments::printSummary() {
  std::string boolStr[2] = {std::string("OFF"), std::string("ON")};
  Logger::info << "General information:" << std::endl;
  Logger::info << "- Output prefix: " << outputPath << std::endl;
  Logger::info << "- Families information: " << familyFilePath << std::endl;
  Logger::info << "- Species tree: " << speciesTree << std::endl;
#ifdef WITH_MPI
  Logger::info << "- MPI Ranks: " << ParallelContext::getSize() << std::endl;
#else
  Logger::info << "- You are running GeneRax without MPI (no parallelization)" << std::endl;
#endif
  Logger::info << "- Random seed: " << seed << std::endl;
  Logger::info << "- Reconciliation model: " << reconciliationModelStr << std::endl;
  if (reconciliationSampleNumber) {
    Logger::info << "- Reconciliation samples: " << reconciliationSampleNumber << std::endl;
  }
  Logger::info << "- DTL rates: "; 
  if (perSpeciesDTLRates) {
    Logger::info << "per species rates" << std::endl;
  } else if (perFamilyDTLRates) {
    Logger::info << "per family rates" << std::endl;
  } else {
    Logger::info << "global rates" << std::endl;
  }
  Logger::info << "- Infer ML reconciliation: " << boolStr[reconcile] << std::endl;
  Logger::info << "- Unrooted reconciliation likelihood: " << boolStr[!rootedGeneTree] << std::endl;
  Logger::info << "- Prune species tree mode: " << boolStr[pruneSpeciesTree] << std::endl;
  Logger::info << std::endl;

  if (speciesStrategy != SpeciesSearchStrategy::SKIP) {
    Logger::info << "Species tree inference information:" << std::endl;
    Logger::info << "- Species tree Strategy: " << ArgumentsHelper::speciesStrategyToStr(speciesStrategy) << std::endl;
    Logger::info << "- Quartet branch supports estimation: " <<  boolStr[quartetSupport] << std::endl;
    Logger::info << "- Branch length estimation" <<  boolStr[estimateSpeciesBranchLenghts] << std::endl;
    Logger::info << "- SPR radius: " <<  speciesSPRRadius << std::endl;
    Logger::info << std::endl;
  } 
    
  if (geneSearchStrategy != GeneSearchStrategy::SKIP) {
    Logger::info << "Gene tree correction information:" << std::endl;  
    Logger::info << "- Gene tree search strategy: " << ArgumentsHelper::strategyToStr(geneSearchStrategy) << std::endl;
    Logger::info << "- Max gene SPR radius: " << maxSPRRadius << std::endl;
    Logger::info << std::endl;
  }
}
