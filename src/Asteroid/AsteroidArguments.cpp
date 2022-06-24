#include "AsteroidArguments.hpp"
#include <IO/Logger.hpp>
   


AsteroidArguments::AsteroidArguments(int iargc, char * iargv[]):
  argc(iargc),
  argv(iargv),
  speciesTreeAlgorithm(SpeciesTreeAlgorithm::User),
  output("Asteroid"),
  minbl(-1.0),
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
    } else if (arg == "-b" || arg == "--min-bl") {
      minbl = atof(argv[++i]);
    } else if (arg == "--seed") {
      seed = atoi(argv[++i]);
    } else if (arg == "--missing-data") {
      missingData = true;
    }
  }
}


void AsteroidArguments::printHelp()
{
  Logger::info << "TODO" << std::endl;
}
