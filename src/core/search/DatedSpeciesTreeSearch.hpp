#pragma once

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class SpeciesSearchState;

class DatedSpeciesTreeSearch {
public:
  static double optimizeDates(SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
      SpeciesSearchState &searchState,
      bool thorough);

};



