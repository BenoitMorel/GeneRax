#pragma once

#include "util/enums.hpp"
#include <string>
#include "IO/ArgumentsHelper.hpp"

  

class GeneTegratorArguments {
public:
  GeneTegratorArguments(int argc, char * argv[]);
  void printHelp();
public:
  int argc;
  char ** argv;
  std::string families;
  std::string speciesTree;
  std::string reconciliationModelStr;
  TransferConstaint transferConstraint;
  OriginationStrategy originationStrategy;
  SpeciesTreeAlgorithm speciesTreeAlgorithm;
  SpeciesSearchStrategy speciesSearchStrategy;
  bool pruneSpeciesTree;
  unsigned int gammaCategories;
  int minCoveredSpecies;
  double trimFamilyRatio;
  unsigned int geneTreeSamples;
  std::string output;
  int seed;
  
  // experimental
  bool randomSpeciesRoot;

};

