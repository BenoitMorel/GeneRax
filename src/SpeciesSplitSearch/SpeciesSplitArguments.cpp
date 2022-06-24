#include "SpeciesSplitArguments.hpp"
#include <IO/Logger.hpp>
   


SpeciesSplitArguments::SpeciesSplitArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  output("SpeciesSplit"),
  seed(123),
  rootedScore(false)
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
    } else if (arg == "-p" || arg == "--prefix") {
      output = std::string(argv[++i]);
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--rooted-score") {
      rootedScore = true;
    }
  }
}


void SpeciesSplitArguments::printHelp()
{
  Logger::info << "TODO" << std::endl;
}
