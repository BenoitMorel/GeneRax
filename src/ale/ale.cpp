#include "AleArguments.hpp"
#include <cstdio>
#include <ccp/ConditionalClades.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "AleOptimizer.hpp"
#include "TrimFamilies.hpp"
#include <util/Paths.hpp>
#include <util/RecModelInfo.hpp>
#include <routines/Routines.hpp>
#include <routines/SlavesMain.hpp>
#include <IO/HighwayCandidateParser.hpp>

void filterInvalidFamilies(Families &families)
{
  Logger::timed << "Filtering families" << std::endl;
  Families validFamilies;
  for (const auto &family: families) {
    std::ifstream is(family.startingGeneTree);
    if (!is || is.peek() == std::ifstream::traits_type::eof()) {
      Logger::error << "Can't open input gene trees for family " << family.name << std::endl;
      continue;
    }
    validFamilies.push_back(family);  
  }
  families = validFamilies;
}

void generateCCPs(const std::string &ccpDir, 
    Families &families, 
    CCPRooting ccpRooting)
{
  ParallelContext::barrier();
  Logger::timed << "Generating ccp files..." << std::endl;
  for (auto &family: families) {
    family.ccp = FileSystem::joinPaths(ccpDir, family.name + ".ccp");
  } 
  auto N = families.size();
  for (auto i = ParallelContext::getBegin(N); i < ParallelContext::getEnd(N); i ++) {
    ConditionalClades ccp(families[i].startingGeneTree, families[i].likelihoodFile, ccpRooting);
    ccp.serialize(families[i].ccp);
  }
  ParallelContext::barrier();
}

void cleanupCCPs(Families &families) 
{
  ParallelContext::barrier();
  Logger::timed << "Cleaning up ccp files..." << std::endl;
  for (const auto &family: families) {
    std::remove(family.ccp.c_str());
  }
  ParallelContext::barrier();
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

void initStartingSpeciesTree(AleArguments &args,
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

void run( AleArguments &args)
{
  Random::setSeed(static_cast<unsigned int>(args.seed));
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  std::string ccpDir = FileSystem::joinPaths(args.output, "ccps");
  FileSystem::mkdir(ccpDir, true);
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "genetegrator"));
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  filterInvalidFamilies(families);
  generateCCPs(ccpDir, families, args.ccpRooting);
  trimFamilies(families, args.minCoveredSpecies, args.trimFamilyRatio);
  if (families.size() == 0) {
    Logger::info << "No valid family, aborting" << std::endl;
    ParallelContext::abort(0);
  }
  initStartingSpeciesTree(args, families);
  
  RecModelInfo info;
  info.pruneSpeciesTree = args.pruneSpeciesTree;
  info.noTL = args.noTL;
  info.model = ArgumentsHelper::strToRecModel(args.reconciliationModelStr); 
  info.transferConstraint = args.transferConstraint;
  info.gammaCategories = args.gammaCategories;
  info.originationStrategy = args.originationStrategy;
  info.fractionMissingFile = args.fractionMissingFile;
  AleOptimizer speciesTreeOptimizer(
      args.speciesTree,
      families,
      info,
      !args.fixRates,
      args.output);
  if (args.randomSpeciesRoot) {
    Logger::timed << "Random root position!" << std::endl;
    speciesTreeOptimizer.randomizeRoot();
  }
  switch (args.speciesSearchStrategy) {
  case SpeciesSearchStrategy::HYBRID:
    speciesTreeOptimizer.optimize();
    break;
  case SpeciesSearchStrategy::EVAL:
    speciesTreeOptimizer.optimizeModelRates(false);
    break;
  case SpeciesSearchStrategy::REROOT:
    speciesTreeOptimizer.optimizeModelRates(true);
    Logger::timed << "First root search, non thorough" << std::endl;
    speciesTreeOptimizer.rootSearch(2, false);
    Logger::timed << "Second root search, thorough" << std::endl;
    speciesTreeOptimizer.rootSearch(2, true);
    break;
  case SpeciesSearchStrategy::SKIP:
    break;
  default:
    assert(false); // not implemented yet
    break;
  }
  if (args.inferSpeciationOrders) {
    speciesTreeOptimizer.optimizeDates();
    speciesTreeOptimizer.getEvaluator().computeLikelihood();
  }
  Logger::timed <<"Sampling reconciled gene trees... (" << args.geneTreeSamples  << " samples)" << std::endl;
  if (args.highways) {
    // let's infer highways of transfers!
    auto highwayOutput = FileSystem::joinPaths(args.output,
      "highway_best_candidates.txt");
    std::vector<ScoredHighway> candidateHighways;
    // initial candidates
    if (args.highwayCandidateFile.size()) {
      // the user sets the candidates
      auto highways = HighwayCandidateParser::parse(args.highwayCandidateFile,
          speciesTreeOptimizer.getSpeciesTree().getTree());
      for (const auto &highway: highways) {
        candidateHighways.push_back(ScoredHighway(highway));
      }
    } else {
      // automatically search for candidates
      speciesTreeOptimizer.getCandidateHighways(candidateHighways, args.highwayCandidatesStep1);
    }
    // first filtering step
    std::vector<ScoredHighway> filteredHighways;
    speciesTreeOptimizer.filterCandidateHighwaysFast(candidateHighways, filteredHighways);
    filteredHighways.resize(std::min(filteredHighways.size(), size_t(args.highwayCandidatesStep2)));
    std::vector<ScoredHighway> bestHighways;
    speciesTreeOptimizer.selectBestHighways(filteredHighways, bestHighways);
    speciesTreeOptimizer.saveBestHighways(bestHighways,
        highwayOutput);
    auto acceptedHighwayOutput = FileSystem::joinPaths(args.output,
      "highway_accepted_highways.txt");
    std::vector<ScoredHighway> acceptedHighways;
    speciesTreeOptimizer.addHighways(bestHighways, acceptedHighways);
    speciesTreeOptimizer.saveBestHighways(acceptedHighways,
        acceptedHighwayOutput);
  }
  speciesTreeOptimizer.optimizeModelRates(true);
  speciesTreeOptimizer.reconcile(args.geneTreeSamples);
  speciesTreeOptimizer.saveSpeciesTree(); 
  if (args.cleanupCCP) {
    cleanupCCPs(families);
  }
  Logger::timed <<"End of the execution" << std::endl;
}

int genetegrator_main(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm); 
  Logger::init();
  Logger::timed << "GeneTegrator v0.0.0" << std::endl; 
  AleArguments args(argc, argv); 
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


