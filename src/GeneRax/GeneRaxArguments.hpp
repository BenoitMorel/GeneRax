#pragma once

#include "enums.hpp"
#include <string>
#include "Arguments.hpp"
using namespace std;

  

class GeneRaxArguments {
public:
   GeneRaxArguments(int argc, char * argv[]);
   void checkInputs();
   void printHelp();
   void printCommand();
   void printSummary();
public:
   int argc;
   char ** argv;
   string families;
   string speciesTree;
   Strategy strategy;
   string reconciliationModelStr;
   RecOpt reconciliationOpt;
   string output;
   bool rootedGeneTree;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
   int maxSPRRadius;
};

