#pragma once

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class AverageStream;

class SpeciesTransferSearch {
public:
  
  static bool transferSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    AverageStream &averageFastDiff,
    double previousBestLL,
    double &newBestLL);


};
