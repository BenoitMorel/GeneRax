#include "SpeciesSPRSearch.hpp"
#include <search/SpeciesSearchCommon.hpp>
#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>

bool SpeciesSPRSearch::SPRRound(SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState,
  unsigned int radius,
  double previousBestLL,
  double &newBestLL)
{
  Logger::timed << "[Species search] Start SPR Round radius=" << radius << std::endl;
  newBestLL = previousBestLL;
  auto hash1 = speciesTree.getNodeIndexHash(); 
  auto supportValues = std::vector<double>();//_getSupport();
  double maxSupport = 0.2; // ignored for now
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(speciesTree, 
      prunes,
      supportValues,
      maxSupport);
  bool better = false;
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      double testedLL = 0.0;
      if (SpeciesSearchCommon::testSPR(speciesTree, evaluation, 
            prune, regraft, searchState, newBestLL, testedLL)) {
        better = true;
        newBestLL = testedLL;
        auto pruneNode = speciesTree.getNode(prune);
        Logger::timed << "\tbetter tree (LL=" 
          << newBestLL << ", hash=" 
          << speciesTree.getHash() 
          << ") "
          << pruneNode->label 
          << " -> " 
          << speciesTree.getNode(regraft)->label
          << std::endl;
        hash1 = speciesTree.getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        if (SpeciesSearchCommon::veryLocalSearch(speciesTree,
            evaluation,
            searchState,
            prune,
            newBestLL,
            testedLL)) {
          newBestLL = testedLL;
        }
      }
    }
  }
  return better;
}

bool SpeciesSPRSearch::SPRSearch(SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState,
  unsigned int radius,
  double previousBestLL,
  double &newBestLL)
{
  
  Logger::info << std::endl;
  Logger::timed << "[Species search]" << " Starting species tree local SPR search, radius=" 
    << radius << ", bestLL=" << previousBestLL  
    << ", hash=" << speciesTree.getHash() << ")" <<std::endl;
  bool better = false;
  while (SPRRound(speciesTree, evaluation, searchState,
        radius, previousBestLL, newBestLL)) {
    previousBestLL = newBestLL;
    better = true;
  }
  Logger::timed << "After normal search: LL=" << newBestLL << std::endl;
  return better;
}



