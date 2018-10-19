
#include "ALEEvaluation.hpp"

ALEEvaluation::ALEEvaluation(const bpp::PhyloTree& speciestree,
  const SpeciesGeneMap& map)
{
  model.setMap(SpeciesGeneMapper::nodeMapsToStrings(map));
  auto species_tree_str = IO::PhyloTreeToNewick(speciestree);
  model.construct_undated(species_tree_str);
}


void ALEEvaluation::setRates(double dupRate, double lossRate)
{
  model.setRates(dupRate, lossRate);
  model.calculate_undatedEs();
}


double ALEEvaluation::evaluate(const bpp::PhyloTree& genetree)
{
  
  model.reset();
  auto gene_tree_str = IO::PhyloTreeToNewick(genetree);
  std::shared_ptr<approx_posterior> ale(new approx_posterior(gene_tree_str));
  std::vector<std::string> gene_tree_strs;
  gene_tree_strs.push_back(gene_tree_str);
  ale->observation(gene_tree_strs);
  return (double)std::log(model.pun(ale, false));
}

