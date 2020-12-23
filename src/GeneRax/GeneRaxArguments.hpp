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
   SpeciesTreeAlgorithm speciesTreeAlgorithm;
   GeneSearchStrategy strategy;
   SpeciesSearchStrategy speciesStrategy;
   std::string reconciliationModelStr;
   std::string output;
   bool perFamilyDTLRates;
   bool rootedGeneTree;
   bool madRooting;
   bool pruneSpeciesTree;
   double supportThreshold;
   unsigned int recRadius;
   bool perSpeciesDTLRates;
   bool useTransferFrequencies;
   bool userDTLRates;
   double dupRate;
   double lossRate;
   double transferRate;
   bool reconcile;
   bool buildSuperMatrix;
   unsigned int reconciliationSamples;
   unsigned int maxSPRRadius;
   double recWeight;
   int seed;
   bool filterFamilies;
   std::string exec;
   bool fractionMissing;


   // species tree search
   bool constrainSpeciesSearch;
   bool rerootSpeciesTree;
   bool estimateSpeciesBranchLenghts;
   unsigned int speciesSPRRadius;
   int speciesInitialFamiliesSubsamples;
   double minGeneBranchLength;
   bool quartetSupport;
   bool quartetSupportAllQuartets;
   int eqpicRadius;
private:

  void init();
};

