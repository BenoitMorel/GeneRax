#include "SpeciesRaxArguments.hpp"
#include <parallelization/ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <limits>
#include <IO/FileSystem.hpp>
#include <sstream>
#include <trees/SpeciesTree.hpp>
#include <trees/PerCoreGeneTrees.hpp>


void initFolders(const std::string &output, std::vector<FamiliesFileParser::FamilyInfo> &families) 
{
  std::string results = FileSystem::joinPaths(output, "results");
  FileSystem::mkdir(results, true);
  for (auto &family: families) {
    FileSystem::mkdir(FileSystem::joinPaths(results, family.name), true);
  }
}


int internal_main(int argc, char** argv, void* comm)
{
  // the order is very important
  ParallelContext::init(comm); 
  Logger::init();
  SpeciesRaxArguments arguments(argc, argv);
  FileSystem::mkdir(arguments.output, true);
  Logger::initFileOutput(FileSystem::joinPaths(arguments.output, "generax"));
  
  arguments.printCommand();
  arguments.printSummary();
  
  RecModel recModel = arguments.reconciliationModel;

  std::vector<FamiliesFileParser::FamilyInfo> initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << std::endl;
  initFolders(arguments.output, initialFamilies);
  
  
  PerCoreGeneTrees geneTrees(initialFamilies); 
  SpeciesTree speciesTree(initialFamilies);
  speciesTree.saveToFile(FileSystem::joinPaths(arguments.output, "starting_species_tree.newick"));
  DTLRates rates(0.1, 0.2, 0.1);
  speciesTree.setRates(rates);
  Logger::info << speciesTree << std::endl;
  //SpeciesTree speciesTree(arguments.speciesTree);
  //SpeciesTreeOperator::changeRoot(speciesTree, 1);
  //SpeciesTreeOperator::changeRoot(speciesTree, 1);
  //Logger::info << "Reconciliation likelihood " << speciesTree.computeReconciliationLikelihood(geneTrees, recModel) << std::endl;
  SpeciesTreeOptimizer::rootExhaustiveSearch(speciesTree, geneTrees, recModel);

  speciesTree.saveToFile(FileSystem::joinPaths(arguments.output, "inferred_species_tree.newick"));
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

