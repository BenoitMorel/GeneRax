#pragma once

#include "util/enums.hpp"
#include <string>
#include "IO/ArgumentsHelper.hpp"

  

class AleArguments {
public:
  AleArguments(int argc, char * argv[]);
  void printHelp();
public:
  
  int argc;
  char ** argv;

  // input data
  std::string families;
  std::string speciesTree;
  std::string reconciliationModelStr;
  
  // model
  TransferConstaint transferConstraint;
  OriginationStrategy originationStrategy;
  bool pruneSpeciesTree;
  unsigned int gammaCategories;
  CCPRooting ccpRooting;

  // search
  SpeciesTreeAlgorithm speciesTreeAlgorithm;
  SpeciesSearchStrategy speciesSearchStrategy;
  bool inferSpeciationOrders;
  bool fixRates;
  bool highways;

  // trimming
  int minCoveredSpecies;
  double trimFamilyRatio;
  
  // output
  unsigned int geneTreeSamples;
  std::string output;
  
  // random seed
  int seed;
  
  // experimental
  bool randomSpeciesRoot;

};
