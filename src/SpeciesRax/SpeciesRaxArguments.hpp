#pragma once

#include "enums.hpp"
#include <string>
using namespace std;

  

class SpeciesRaxArguments {
public:
   SpeciesRaxArguments(int argc, char * argv[]);
   void checkInputs();
   void printHelp();
   void printCommand();
   void printSummary();
public:
   int argc;
   char ** argv;
   string families;
   string speciesTree;
   RecModel reconciliationModel;
   RecOpt reconciliationOpt;
   string output;
   bool rootedGeneTree;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
};

