#pragma once

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

  double computeLogLikelihood();

private:
  PLLRootedTree _speciesTree;
  PerCoreMultiEvaluation _evaluations;

};
