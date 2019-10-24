#include "SlavesMain.hpp"
#include <IO/Logger.hpp>

#include <string>
#include <cassert>

#include <routines/RaxmlSlave.hpp>
#include <routines/GeneRaxSlave.hpp>

bool SlavesMain::isSlave(int argc, char** argv)
{
  if (argc < 2) {
    return false;
  }
  std::string key(argv[1]);
  return key == "optimizeGeneTrees" || key == "raxmlLight";
}

extern "C" int static_scheduled_main(int argc, char** argv, void* comm)
{
  Logger::enableLogFile(false);
  std::string key(argv[1]);
  int res = 1;
  if (key == "optimizeGeneTrees") {
    res = GeneRaxSlave::optimizeGeneTreesMain(argc, argv, comm);   
  } else if (key == "raxmlLight") {
    res = RaxmlSlave::runRaxmlOptimization(argc, argv, comm);
  } else {
    assert(0);
  }
  Logger::enableLogFile(true);
  return res;
}



