#pragma once

#include "util/enums.hpp"
#include <string>
#include "IO/Arguments.hpp"


  

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
   std::string execPath;
   std::string families;
   std::string speciesTree;
   Strategy strategy;
   std::string reconciliationModelStr;
   RecOpt reconciliationOpt;
   std::string output;
   bool perFamilyDTLRates;
   unsigned int duplicates;
   unsigned int initStrategies;
   bool rootedGeneTree;
   bool pruneSpeciesTree;
   unsigned int recRadius;
   bool perSpeciesDTLRates;
   bool useTransferFrequencies;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
   int maxSPRRadius;
   double recWeight;
   unsigned int seed;
   std::string exec;
   
   // species tree search
   bool optimizeSpeciesTree;
   unsigned int speciesFastRadius;
   unsigned int speciesSlowRadius;

};

