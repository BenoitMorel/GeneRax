#pragma once

#include <trees/DatedTree.hpp>

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class SpeciesSearchState;

struct ScoredBackup {
  DatedTree::Backup backup;
  double score;

  ScoredBackup(): score(0.0) {}

  ScoredBackup(DatedTree &datedTree, double score):
    backup(datedTree.getBackup()), score(score) 
  {}

  bool operator < (const ScoredBackup &other) {
    return score < other.score;
  }
};
using ScoredBackups = std::vector<ScoredBackup>;

class DatedSpeciesTreeSearch {
public:
  /**
   * Searchf or the speciation order (dating) that optimizes
   * the score returned by evaluation. If the score gets higher than
   * searchState.bestLL and if searchState.pathToBestSpeciesTree
   * is set, saves the new best tree and update bestLL
   *
   * If thorough is not set, we only apply one naive round.
   * Otherwise, we conduct a more thorough search
   */
  static double optimizeDates(SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
      SpeciesSearchState &searchState,
      double currentLL,
      bool thorough);

  /**
   *
   */
  static ScoredBackups optimizeDatesFromReconciliation(SpeciesTree &speciesTree,
      SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
      unsigned int searches);

};



