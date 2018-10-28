
#include "ALEEvaluation.hpp"
#include <cmath>

ALEEvaluation::ALEEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& map):
  firstCall(true)
{
  undatedDLModel.setGeneSpeciesMap(map);
  undatedDLModel.setSpeciesTree(speciesTree);
}


void ALEEvaluation::setRates(double dupRate, double lossRate)
{
  undatedDLModel.setRates(dupRate, lossRate);
}


double ALEEvaluation::evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  if (firstCall) {
    undatedDLModel.setInitialGeneTree(treeinfo);
  }
  return (double)std::log(undatedDLModel.pun(treeinfo));
}

