#pragma once

#include <vector>
#include <likelihoods/LibpllEvaluation.hpp>

class PLLRootedTree;
class PLLUnrootedTree;

class PolytomySolver {
public:  
  static void solveSimpleInterface(
      PLLRootedTree &speciesTree,
      std::map<std::string, unsigned int> &speciesLabelsToSolve
      );


};
