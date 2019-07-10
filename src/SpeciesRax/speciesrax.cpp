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
  
  std::vector<FamiliesFileParser::FamilyInfo> initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << std::endl;
  initFolders(arguments.output, initialFamilies);
  
  
  PerCoreGeneTrees geneTrees(initialFamilies); 
  SpeciesTree speciesTree(arguments.speciesTree);
  DTLRates rates(0.1, 0.2, 0.1);
  speciesTree.setRates(rates);
  speciesTree.changeRoot(true, false);
  Logger::info << "Tree: " << std::endl << speciesTree << std::endl;
  Logger::info << "Reconciliation likelihood " << speciesTree.computeReconciliationLikelihood(geneTrees, UndatedDTL) << std::endl;
  
  
  
  
  bool left1 = false;
  bool left2 = true;
  if (speciesTree.canChangeRoot(left1, left2)) {
    speciesTree.changeRoot(left1, left2);
    //Logger::info << "Tree: " << std::endl << speciesTree << std::endl;
    Logger::info << "Reconciliation likelihood " << speciesTree.computeReconciliationLikelihood(geneTrees, UndatedDTL) << std::endl;
  } else {
    Logger::info << "Cannot change root" << std::endl;
  }
  left1 = !left1;
  left2 = !left2;
  if (speciesTree.canChangeRoot(left1, left2)) {
    speciesTree.changeRoot(left1, left2);
    Logger::info << "Tree: " << std::endl << speciesTree << std::endl;
    Logger::info << "Reconciliation likelihood " << speciesTree.computeReconciliationLikelihood(geneTrees, UndatedDTL) << std::endl;
  } else {
    Logger::info << "Cannot change root" << std::endl;
  }
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

