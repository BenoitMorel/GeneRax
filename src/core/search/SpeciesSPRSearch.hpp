#pragma once

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class AverageStream;


class SpeciesSPRSearch {
public:
  static bool SPRRound(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    AverageStream &averageFastDiff,
    unsigned int radius,
    double previousBestLL,
    double &newBestLL);
  
  static bool SPRSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    AverageStream &averageFastDiff,
    unsigned int radius,
    double previousBestLL,
    double &newBestLL);

};

