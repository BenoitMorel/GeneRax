#include "GeneTegratorArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "GTSpeciesTreeOptimizer.hpp"
#include "TrimFamilies.hpp"
#include <util/Paths.hpp>
#include <util/RecModelInfo.hpp>
#include <routines/Routines.hpp>
#include <routines/SlavesMain.hpp>


void filterInvalidFamilies(Families &families)
{
  Logger::timed << "Filtering families" << std::endl;
  Families validFamilies;
  for (const auto &family: families) {
    std::ifstream is(family.startingGeneTree);
    if (!is || is.peek() == std::ifstream::traits_type::eof()) {
      Logger::info << "Invalid family " << family.name << std::endl;
      continue;
    }
    validFamilies.push_back(family);  
  }
  families = validFamilies;
}

void trimFamilies(Families &families, int minSpecies, double trimRatio) 
{
  Logger::timed << "Families: " << families.size() << std::endl;
  if (minSpecies != -1) {
    Logger::timed << "Triming families covering less than " << minSpecies << " species " << std::endl;
    TrimFamilies::trimMinSpeciesCoverage(families, minSpecies);
    Logger::timed << "Families: " << families.size() << std::endl;
  }
  if (trimRatio < 1.0) {
    Logger::timed << "Trimming families with too many clades (keeping " 
      << trimRatio * 100.0 << "\% of the families) " << std::endl;
    TrimFamilies::trimHighCladesNumber(families, trimRatio);
  }
  Logger::timed << "Families: " << families.size() << std::endl;
}

void initStartingSpeciesTree(GeneTegratorArguments &args,
    Families &families)
{
  Logger::timed << "Initializing starting species tree..." << std::endl;
  auto startingSpeciesTree = Paths::getSpeciesTreeFile(
      args.output, 
      "starting_species_tree.newick");
  std::unique_ptr<PLLRootedTree> speciesTree(nullptr);
  if (args.speciesTreeAlgorithm == SpeciesTreeAlgorithm::User) {
    unsigned int canRead = 1;
    if (ParallelContext::getRank() == 0) {
      try {
        SpeciesTree reader(args.speciesTree);
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
    PLLRootedTree::labelRootedTree(args.speciesTree, startingSpeciesTree);
  } else {
    Routines::computeInitialSpeciesTree(families,
        args.output,
        args.speciesTreeAlgorithm)->save(startingSpeciesTree);

  }
  ParallelContext::barrier();
  args.speciesTree = Paths::getSpeciesTreeFile(args.output, "inferred_species_tree.newick");
  if (ParallelContext::getRank() == 0) {
    SpeciesTree copy(startingSpeciesTree); 
    copy.getTree().save(args.speciesTree);
  }
  ParallelContext::barrier();
  Logger::timed << "Finished starting species tree initialization" << std::endl;
}

void run( GeneTegratorArguments &args)
{
  Random::setSeed(static_cast<unsigned int>(args.seed));
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "genetegrator"));
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  filterInvalidFamilies(families);
  trimFamilies(families, args.minCoveredSpecies, args.trimFamilyRatio);
  if (families.size() == 0) {
    Logger::info << "No valid family, aborting" << std::endl;
    ParallelContext::abort(0);
  }
  initStartingSpeciesTree(args, families);
  
  RecModelInfo info;
  info.pruneSpeciesTree = args.pruneSpeciesTree;
  info.model = ArgumentsHelper::strToRecModel(args.reconciliationModelStr); 
  info.transferConstraint = args.transferConstraint;
  info.gammaCategories = args.gammaCategories;
  info.originationStrategy = args.originationStrategy;
  GTSpeciesTreeOptimizer speciesTreeOptimizer(
      args.speciesTree,
      families,
      info,
      args.output);
  if (args.randomSpeciesRoot) {
    Logger::timed << "Random root position!" << std::endl;
    speciesTreeOptimizer.randomizeRoot();
    Logger::timed << "New ll=" << speciesTreeOptimizer.getEvaluator().computeLikelihood() << std::endl;
  }
  switch (args.speciesSearchStrategy) {
  case SpeciesSearchStrategy::HYBRID:
    speciesTreeOptimizer.optimize();
    break;
  case SpeciesSearchStrategy::EVAL:
    
    //speciesTreeOptimizer.optimizeModelRates(false);
    //speciesTreeOptimizer.optimizeDates();
    //speciesTreeOptimizer.optimizeModelRates(false);
    speciesTreeOptimizer.optimizeModelRates(true);
    //speciesTreeOptimizer.optimizeDates();
    //Logger::timed << "First root search, non thorough" << std::endl;
    //speciesTreeOptimizer.rootSearch(10, false);
    //Logger::timed << "Second root search, thorough" << std::endl;
    //speciesTreeOptimizer.rootSearch(2, true);
    break;
  case SpeciesSearchStrategy::SKIP:
    break;
  default:
    assert(false); // not implemented yet
    break;
  }
  Logger::timed <<"Sampling reconciled gene trees... (" << args.geneTreeSamples  << " samples)" << std::endl;
  speciesTreeOptimizer.reconcile(args.geneTreeSamples);
  speciesTreeOptimizer.saveSpeciesTree(); 
  Logger::timed <<"End of the execution" << std::endl;
}

int genetegrator_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm); 
  Logger::init();
  Logger::timed << "GeneTegrator v0.0.0" << std::endl; 
  GeneTegratorArguments args(argc, argv); 
  run(args);
  Logger::close();
  ParallelContext::finalize();
  return 0;
}

int internal_main(int argc, char** argv, void* comm)
{
  if (SlavesMain::isSlave(argc, argv)) {
    int slaveComm = -1; 
    return static_scheduled_main(argc, argv, &slaveComm);
  } else {
    return genetegrator_main(argc, argv, comm);
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


