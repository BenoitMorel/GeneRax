#include "GeneTegratorArguments.hpp"
#include <IO/Logger.hpp>
   


GeneTegratorArguments::GeneTegratorArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  pruneSpeciesTree(false),
  output("GeneTegrator"),
  seed(123)
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
    } else if (arg == "--prune-species-tree") {
      pruneSpeciesTree = true;
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    }
  }
}


void GeneTegratorArguments::printHelp()
{
  Logger::info << "TODO" << std::endl;
}
