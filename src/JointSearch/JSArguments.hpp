#pragma once

#include "enums.hpp"
#include <string>
using namespace std;

  

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
   string geneTree;
   string alignment;
   string speciesTree;
   string geneSpeciesMap;
   Strategy strategy;
   RecModel reconciliationModel;
   RecOpt reconciliationOpt;
   string libpllModel;
   string output;
   bool check;
   bool rootedGeneTree;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
};

