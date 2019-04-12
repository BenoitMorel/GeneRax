#include "SpeciesRaxArguments.hpp"
#include <ParallelContext.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <limits>
#include <IO/FileSystem.hpp>
#include <sstream>

using namespace std;

void initFolders(const string &output, vector<FamiliesFileParser::FamilyInfo> &families) 
{
  string results = FileSystem::joinPaths(output, "results");
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

  vector<FamiliesFileParser::FamilyInfo> initialFamilies = FamiliesFileParser::parseFamiliesFile(arguments.families);
  Logger::info << "Number of gene families: " << initialFamilies.size() << endl;
  initFolders(arguments.output, initialFamilies);
  ParallelContext::finalize();
  return 0;
}



int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

