
#include "ALEEvaluation.hpp"

ALEEvaluation::ALEEvaluation(const bpp::PhyloTree& speciestree,
    pll_rtree_t *pllSpeciesTree,
  const SpeciesGeneMap& map)
{
  undatedDLModel.setMap(SpeciesGeneMapper::nodeMapsToStrings(map));
  auto species_tree_str = IO::PhyloTreeToNewick(speciestree);
  undatedDLModel.setSpeciesTree(species_tree_str, pllSpeciesTree);
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

