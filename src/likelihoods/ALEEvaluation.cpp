
#include "ALEEvaluation.hpp"
#include <cmath>

ALEEvaluation::ALEEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& map,
  bool transfers):
  transfers(transfers),
  firstCall(true)
{
  if (transfers) {
    undatedDTLModel.setGeneSpeciesMap(map);
    undatedDTLModel.setSpeciesTree(speciesTree);
  } else {
    undatedDLModel.setGeneSpeciesMap(map);
    undatedDLModel.setSpeciesTree(speciesTree);
  }
}


void ALEEvaluation::setRates(double dupRate, double lossRate,
  double transferRate)
{
  if (transfers) {
    undatedDTLModel.setRates(dupRate, lossRate, transferRate);
  } else {
    assert(transferRate == 0.0);
    undatedDLModel.setRates(dupRate, lossRate);
  }
}


double ALEEvaluation::evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  if (transfers) {
    if (firstCall) {
      undatedDTLModel.setInitialGeneTree(treeinfo);
    }
    return (double)std::log(undatedDTLModel.pun(treeinfo));
  } else {
    if (firstCall) {
      undatedDLModel.setInitialGeneTree(treeinfo);
    }
    return (double)std::log(undatedDLModel.pun(treeinfo));
  }
}

