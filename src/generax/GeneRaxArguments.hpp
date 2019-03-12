#pragma once

#include <string>
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
   string strategy;
   string reconciliationModel;
   string reconciliationOpt;
   string libpllModel;
   string output;
   bool check;
   bool rootedGeneTree;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
};

