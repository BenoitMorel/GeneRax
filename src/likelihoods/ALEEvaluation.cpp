
#include "ALEEvaluation.hpp"

ALEEvaluation::ALEEvaluation(const bpp::PhyloTree& speciestree,
  const SpeciesGeneMap& map)
{
  model.setMap(SpeciesGeneMapper::nodeMapsToStrings(map));
  auto species_tree_str = IO::PhyloTreeToNewick(speciestree);
  model.construct_undated(species_tree_str, "");
}


void ALEEvaluation::setRates(double dupRate, double lossRates)
{

}


  double ALEEvaluation::evaluate(const bpp::PhyloTree& genetree,
    const bpp::PhyloTree& speciestree,
    const SpeciesGeneMap& map,
    long double beta,
    long double O_R,
    long double delta,
    long double tau,
    long double lambda
    ) {
  model.reset();
  model.set_model_parameter("seq_beta", beta);
  model.set_model_parameter("O_R", O_R);
  model.set_model_parameter("delta", delta);
  model.set_model_parameter("tau", tau);
  model.set_model_parameter("lambda", lambda);
  
  
  model.calculate_undatedEs();


  auto gene_tree_str = IO::PhyloTreeToNewick(genetree);
  std::shared_ptr<approx_posterior> ale(new approx_posterior(gene_tree_str));
  std::vector<std::string> gene_tree_strs;
  gene_tree_strs.push_back(gene_tree_str);
  ale->observation(gene_tree_strs);

  double loglK = (double)std::log(model.pun(ale, false));

  return loglK;
}

