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
   std::string familyFilePath;
   std::string speciesTree;
   GeneSearchStrategy geneSearchStrategy;
   std::string reconciliationModelStr;
   std::string outputPath;
   bool perFamilyDTLRates;
   bool rootedGeneTree;
   bool madRooting;
   bool pruneSpeciesTree;
   double supportThreshold;
   unsigned int recRadius;
   bool perSpeciesDTLRates;
   bool userDTLRates;
   TransferConstaint transferConstraint;
   bool noDup;
   double dupRate;
   double lossRate;
   double transferRate;
   size_t gammaCategories;
   OriginationStrategy originationStrategy;
   bool reconcile;
   bool buildSuperMatrix;
   unsigned int reconciliationSampleNumber;
   unsigned int maxSPRRadius;
   double recWeight;
   int seed;
   bool filterFamilies;
   std::string exec;
   std::string fractionMissingFile;


   // species tree search
   SpeciesTreeAlgorithm speciesTreeAlgorithm;
   SpeciesSearchStrategy speciesStrategy;
   bool constrainSpeciesSearch;
   bool rerootSpeciesTree;
   bool estimateSpeciesBranchLenghts;
   unsigned int speciesSPRRadius;
   unsigned int speciesSmallRootRadius;
   unsigned int speciesBigRootRadius;
   double minGeneBranchLength;
   bool quartetSupport;
   bool quartetSupportAllQuartets;
   int eqpicRadius;
   
   // hacky stuff
   bool generateFakeAlignments;
private:

  void init();
};

