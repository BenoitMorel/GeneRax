#pragma once

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class AverageStream;
class SpeciesSearchState;

class SpeciesTransferSearch {
public:
  
  static bool transferSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState);


};
