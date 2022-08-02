#pragma once

#include <vector>
#include <likelihoods/LibpllEvaluation.hpp>
#include <map>

class PLLRootedTree;
class PLLUnrootedTree;

class PolytomySolver {
public:  
  static void solveSimpleInterface(
      PLLRootedTree &speciesTree,
      std::map<std::string, unsigned int> &speciesLabelsToSolve
      );


};
