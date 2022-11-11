#pragma once

#include <string>

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class AverageStream;
class SpeciesSearchState;

class SpeciesTransferSearch {
public:
  
  static bool transferSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    const std::string &outputDir);


};
