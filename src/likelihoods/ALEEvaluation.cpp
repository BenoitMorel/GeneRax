
#include "ALEEvaluation.hpp"

ALEEvaluation::ALEEvaluation(pll_rtree_t *speciesTree,
  const SpeciesGeneMap& map)
{
  undatedDLModel.setMap(SpeciesGeneMapper::nodeMapsToStrings(map));
  undatedDLModel.setSpeciesTree(speciesTree);
}


void ALEEvaluation::setRates(double dupRate, double lossRate)
{
  undatedDLModel.setRates(dupRate, lossRate);
}


double ALEEvaluation::evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  undatedDLModel.reset();
  return (double)std::log(undatedDLModel.pun(treeinfo));
}

