#pragma once

#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "UndatedDLMultiModel.hpp"
#include <trees/PLLRootedTree.hpp>
#include <memory>
#include <vector>

using MultiEvaluation = UndatedDLMultiModel<ScaledValue>;
using MultiEvaluationPtr = 
  std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluation = std::vector<MultiEvaluationPtr>;

class GTSpeciesTreeOptimizer {
public:
  GTSpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &families, 
      const std::string &outputDir);

  double computeRecLikelihood();
  double sprSearch(unsigned int radius);

private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  PerCoreMultiEvaluation _evaluations;
  std::string _outputDir;
  double _bestRecLL;

  double fastSPRRound(unsigned int radius);
  bool testPruning(unsigned int prune,
    unsigned int regraft);
  void newBestTreeCallback(double newLL);
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};
