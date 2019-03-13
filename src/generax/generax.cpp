#include "GeneRaxArguments.hpp"
#include <ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <limits>
#include <PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
using namespace std;


void optimizeGeneTrees(vector<FamiliesFileParser::FamilyInfo> &families) 
{

}

int internal_main(int argc, char** argv, void* comm)
{
  // the order is very important
  ParallelContext::init(comm); 
  Logger::init();
  GeneRaxArguments arguments(argc, argv);
  Logger::initFileOutput(arguments.output);
  
  arguments.printCommand();
  arguments.printSummary();

  vector<FamiliesFileParser::FamilyInfo> families = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << families.size() << endl;

  PerCoreGeneTrees geneTrees(families);
  pll_rtree_t *speciesTree = LibpllEvaluation::readRootedFromFile(arguments.speciesTree); 
  DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, arguments.reconciliationModel);
  Logger::timed << "End of GeneRax execution" << endl;
  ParallelContext::finalize();
  return 0;
}


#ifdef JOINTSEARCH_BUILD_AS_LIB

extern "C" int dll_main(int argc, char** argv, void* comm)
{
  return internal_main(argc, argv, comm);
}

#else

int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

#endif

