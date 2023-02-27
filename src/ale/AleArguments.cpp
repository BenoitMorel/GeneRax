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
  gammaCategories(1),
  ccpRooting(CCPRooting::UNIFORM),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  speciesSearchStrategy(SpeciesSearchStrategy::HYBRID),
  inferSpeciationOrders(false),
  fixRates(false),
  highways(false),
  minCoveredSpecies(4),
  trimFamilyRatio(1.0),
  geneTreeSamples(0),
  output("GeneTegrator"),
  seed(123),
  randomSpeciesRoot(false)
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
    } else if (arg == "--fix-rates") {
      fixRates = true;
    } else if (arg == "--highways") {
      highways = true;
    } else if (arg == "--transfer-constraint") {
      transferConstraint = ArgumentsHelper::strToTransferConstraint(std::string(argv[++i]));
    } else if (arg == "--origination") {
      originationStrategy = Enums::strToOrigination(std::string(argv[++i]));
    } else if (arg == "--prune-species-tree") {
      pruneSpeciesTree = true;
    } else if (arg == "--trim-ratio") {
      trimFamilyRatio = atof(argv[++i]);
    } else if (arg == "--min-covered-species") {
      minCoveredSpecies = atof(argv[++i]);
    } else if (arg == "--gamma-categories") {
      gammaCategories = atoi(argv[++i]);
    } else if (arg == "--gene-tree-rooting") {
      ccpRooting = ArgumentsHelper::strToCCPRooting(std::string(argv[++i]));
    } else if (arg == "--gene-tree-samples") {
      geneTreeSamples = atoi(argv[++i]);
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--random-species-root") {
      randomSpeciesRoot = true;
    } else {
      std::cerr << "Unknown argument " << arg << std::endl;
    }
  }
}


void AleArguments::printHelp()
{
  Logger::info << "TODO" << std::endl;
}
