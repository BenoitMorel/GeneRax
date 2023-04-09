#include "AleArguments.hpp"
#include <IO/Logger.hpp>
#include <climits>  


AleArguments::AleArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  reconciliationModelStr("UndatedDTL"),
  transferConstraint(TransferConstaint::PARENTS),
  originationStrategy(OriginationStrategy::UNIFORM),
  pruneSpeciesTree(false),
  noTL(false),
  gammaCategories(1),
  ccpRooting(CCPRooting::UNIFORM),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  speciesSearchStrategy(SpeciesSearchStrategy::HYBRID),
  inferSpeciationOrders(false),
  fixRates(false),
  skipThoroughRates(false),
  highways(false),
  highwayCandidatesStep1(100),
  highwayCandidatesStep2(25),
  minCoveredSpecies(4),
  trimFamilyRatio(1.0),
  geneTreeSamples(0),
  output("GeneTegrator"),
  cleanupCCP(true),
  seed(123),
  randomSpeciesRoot(false),
  verboseOptRates(false)
{
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      printHelp();
      ParallelContext::abort(0);
    } else if (arg == "-f" || arg == "--families") {
      families = std::string(argv[++i]);
    } else if (arg == "-s" || arg == "--species-tree") {
      speciesTree = std::string(argv[++i]);
      speciesTreeAlgorithm = Enums::strToSpeciesTree(speciesTree);
    } else if (arg == "--species-search") {
      speciesSearchStrategy = ArgumentsHelper::strToSpeciesSearchStrategy(std::string(argv[++i]));
    } else if (arg == "--infer-speciation-orders") {
      inferSpeciationOrders = true;
    } else if (arg == "-r" || arg == "--rec-model") {
      reconciliationModelStr = std::string(argv[++i]);
    } else if (arg == "--skip-thorough-rates") {
      skipThoroughRates = true;
    } else if (arg == "--fix-rates") {
      fixRates = true;
    } else if (arg == "--highways") {
      highways = true;
    } else if (arg == "--highway-candidates-file") {
      highwayCandidateFile= std::string(argv[++i]);
    } else if (arg == "--highway-candidates-step1") {
      highwayCandidatesStep1 = atoi(argv[++i]);
    } else if (arg == "--highway-candidates-step2") {
      highwayCandidatesStep2 = atoi(argv[++i]);
    } else if (arg == "--transfer-constraint") {
      transferConstraint = ArgumentsHelper::strToTransferConstraint(std::string(argv[++i]));
    } else if (arg == "--origination") {
      originationStrategy = Enums::strToOrigination(std::string(argv[++i]));
    } else if (arg == "--prune-species-tree") {
      pruneSpeciesTree = true;
    } else if (arg == "--no-tl") {
      noTL = true;
    } else if (arg == "--trim-ratio") {
      trimFamilyRatio = atof(argv[++i]);
    } else if (arg == "--min-covered-species") {
      minCoveredSpecies = atof(argv[++i]);
    } else if (arg == "--gamma-categories") {
      gammaCategories = atoi(argv[++i]);
    } else if (arg == "--gene-tree-rooting") {
      ccpRooting = ArgumentsHelper::strToCCPRooting(std::string(argv[++i]));
    } else if (arg == "--fraction-missing") {
      fractionMissingFile = argv[++i];
    } else if (arg == "--species-categories") {
      speciesCategoryFile = argv[++i];
    } else if (arg == "--gene-tree-samples") {
      geneTreeSamples = atoi(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--skip-cleanup-ccp") {
      cleanupCCP = false;
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--random-species-root") {
      randomSpeciesRoot = true;
    } else if (arg == "--verbose-opt-rates") {
      verboseOptRates = true;
    } else {
      std::cerr << "Unknown argument " << arg << std::endl;
    }
  }
}


void AleArguments::printHelp()
{
  Logger::info << "TODO" << std::endl;
}
