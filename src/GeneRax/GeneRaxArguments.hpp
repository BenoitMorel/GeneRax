#pragma once

#include "util/enums.hpp"
#include <string>
#include "IO/ArgumentsHelper.hpp"


  

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
   SpeciesStrategy speciesStrategy;
   std::string reconciliationModelStr;
   RecOpt reconciliationOpt;
   std::string output;
   bool perFamilyDTLRates;
   bool rootedGeneTree;
   bool pruneSpeciesTree;
   double supportThreshold;
   unsigned int recRadius;
   bool perSpeciesDTLRates;
   bool useTransferFrequencies;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
   bool optimizeGeneTrees;
   bool reconcile;
   unsigned int reconciliationSamples;
   unsigned int maxSPRRadius;
   double recWeight;
   int seed;
   std::string exec;
   
   // species tree search
   bool optimizeSpeciesTree;
   unsigned int speciesFastRadius;
   unsigned int speciesSlowRadius;
private:

  void init();
};

