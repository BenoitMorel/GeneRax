#include "MiniBMEArguments.hpp"
#include <IO/Logger.hpp>
   


MiniBMEArguments::MiniBMEArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  output("MiniBME"),
  seed(123),
  missingData(false)
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
    } else if (arg == "--missing-data") {
      missingData = true;
    }
  }
}


void MiniBMEArguments::printHelp()
{
  Logger::info << "TODO" << std::endl;
}
