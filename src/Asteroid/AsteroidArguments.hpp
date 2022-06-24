#pragma once

#include "util/enums.hpp"
#include <string>
#include "IO/ArgumentsHelper.hpp"


  

class AsteroidArguments {
public:
  AsteroidArguments(int argc, char * argv[]);
  void printHelp();
public:
   int argc;
   char ** argv;
   std::string families;
   std::string speciesTree;
   SpeciesTreeAlgorithm speciesTreeAlgorithm;
   std::string output;
   double minbl;
   int seed;
   bool missingData;

};

