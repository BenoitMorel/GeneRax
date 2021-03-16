#include "GeneTegratorArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "GTSpeciesTreeOptimizer.hpp"
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
    LibpllParsers::labelRootedTree(args.speciesTree, startingSpeciesTree);
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

}

void run( GeneTegratorArguments &args)
{
  Logger::info << "Mkdir " << args.output << std::endl;
  FileSystem::mkdir(args.output, true);
  FileSystem::mkdir(args.output + "/species_trees", true);
  Logger::initFileOutput(FileSystem::joinPaths(args.output, "genetegrator"));
  auto families = FamiliesFileParser::parseFamiliesFile(args.families);
  filterInvalidFamilies(families);
  initStartingSpeciesTree(args, families);
  
  RecModelInfo info;
  info.pruneSpeciesTree = args.pruneSpeciesTree;
  info.model = ArgumentsHelper::strToRecModel(args.reconciliationModelStr); 
  GTSpeciesTreeOptimizer speciesTreeOptimizer(
      args.speciesTree,
      families,
      info,
      args.output);
  auto ll = speciesTreeOptimizer.computeRecLikelihood();
  Logger::info << "total ll=" << ll << std::endl;
  ll = speciesTreeOptimizer.sprSearch(1);
  
  Logger::info << "total ll=" << ll << std::endl;
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


