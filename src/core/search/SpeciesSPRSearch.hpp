#pragma once

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class AverageStream;
class SpeciesSearchState;


class SpeciesSPRSearch {
public:
  static bool SPRRound(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int radius);
  
  static bool SPRSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int radius);

};

