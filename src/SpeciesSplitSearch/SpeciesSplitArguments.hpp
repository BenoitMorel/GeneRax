#pragma once

#include "util/enums.hpp"
#include <string>
#include "IO/ArgumentsHelper.hpp"


  

class SpeciesSplitArguments {
public:
  SpeciesSplitArguments(int argc, char * argv[]);
  void printHelp();
public:
   int argc;
   char ** argv;
   std::string families;
   std::string speciesTree;
   SpeciesTreeAlgorithm speciesTreeAlgorithm;
   std::string output;
   int seed;
   bool rootedScore;

};

