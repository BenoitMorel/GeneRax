#pragma once

#include "util/enums.hpp"
#include <string>


  

class JSArguments {
public:
   JSArguments(int argc, char * argv[]);
   void checkInputs();
   void printHelp();
   void printCommand();
   void printSummary();
public:
   int argc;
   char ** argv;
   std::string geneTree;
   std::string alignment;
   std::string speciesTree;
   std::string geneSpeciesMap;
   GeneSearchStrategy strategy;
   RecModel reconciliationModel;
   RecOpt reconciliationOpt;
   std::string libpllModel;
   std::string output;
   bool check;
   bool rootedGeneTree;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
};

