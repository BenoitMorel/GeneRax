#pragma once

#include <string>
using namespace std;

  

class GecoArguments {
public:
   GecoArguments(int argc, char * argv[]);
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

