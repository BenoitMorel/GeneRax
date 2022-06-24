#include "SpeciesSPRSearch.hpp"
#include <search/SpeciesSearchCommon.hpp>
#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>

bool SpeciesSPRSearch::SPRRound(SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState,
  unsigned int radius)
{
  Logger::timed << "[Species search] Start SPR Round radius=" << radius << std::endl;
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
      if (SpeciesSearchCommon::testSPR(speciesTree, evaluation, 
            searchState, prune, regraft)) {
        better = true;
        auto pruneNode = speciesTree.getNode(prune);
        Logger::timed << "\tbetter tree (LL=" 
          << searchState.bestLL << ", hash=" 
          << speciesTree.getHash() 
          << ") "
          << pruneNode->label 
          << " -> " 
          << speciesTree.getNode(regraft)->label
          << std::endl;
        hash1 = speciesTree.getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        SpeciesSearchCommon::veryLocalSearch(speciesTree,
            evaluation,
            searchState,
            prune);
      }
    }
  }
  return better;
}

bool SpeciesSPRSearch::SPRSearch(SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState,
  unsigned int radius)
{
  
  Logger::info << std::endl;
  Logger::timed << "[Species search]" << " Starting species tree local SPR search, radius=" 
    << radius << ", bestLL=" << searchState.bestLL 
    << ", hash=" << speciesTree.getHash() << ")" <<std::endl;
  bool better = false;
  while (SPRRound(speciesTree, evaluation, searchState,
        radius)) {
    better = true;
  }
  Logger::timed << "After normal search: LL=" << searchState.bestLL << std::endl;
  return better;
}



