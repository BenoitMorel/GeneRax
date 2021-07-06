#pragma once

#include <search/SpeciesSearchCommon.hpp>

class SpeciesSearchState;

class SpeciesRootSearch {
public:
  /**
   *  Search for the ML root for the current
   *  species tree topology
   *  
   *  rootLikelihoods and treePerFamLLVec 
   *  will only be filled if not NULL
   */
  static double rootSearch(
      SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
      SpeciesSearchState &searchState,
      unsigned int maxDepth,
      RootLikelihoods *rootLikelihoods = nullptr,
      TreePerFamLLVec *treePerFamLLVec = nullptr); 
};

