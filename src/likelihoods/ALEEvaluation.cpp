
#include "ALEEvaluation.hpp"
#include <cmath>

#include "likelihoods/dated_ale/ALE.h"
#include "likelihoods/dated_ale/exODT.h"


ALEEvaluation::ALEEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& map):
  firstCall(true)
{
#ifdef UNDATED
  undatedDLModel.setGeneSpeciesMap(map);
  undatedDLModel.setSpeciesTree(speciesTree);
#else
  speciesTreeStr_ = string(pll_rtree_export_newick(speciesTree->root, 0));
   
#endif
}


void ALEEvaluation::setRates(double dupRate, double lossRate)
{
#ifdef UNDATED
  undatedDLModel.setRates(dupRate, lossRate);
#else
  dupRate_ = dupRate;
  lossRate_ = lossRate;
  transferRate_ = 0;
#endif
}



approx_posterior * observe_ALE_from_string(string tree)
{
  vector<string> trees;
  trees.push_back(tree);
  approx_posterior* ale=new approx_posterior(trees[0]);// NO del-loc
  ale->observation(trees,false);
  return ale;
}

double evalDated(const string &geneTreeStr,
  const string &speciesTreeStr,
  double dupRate,
  double lossRate,
  double transferRate)
{/*
  cout << "eval dated " << endl;
  cout << geneTreeStr << endl;
  cout << speciesTreeStr << endl;
  */
  cout << dupRate << " " << lossRate << " " << transferRate << endl;
  auto ale = observe_ALE_from_string(geneTreeStr);
  auto model = new exODT_model();
  model->set_model_parameter("min_D",3);
  model->set_model_parameter("grid_delta_t",0.05);
  model->construct(speciesTreeStr);
  
  model->set_model_parameter("event_node",0);
  model->set_model_parameter("leaf_events",1);

  model->set_model_parameter("delta", dupRate);
  model->set_model_parameter("tau", transferRate);
  model->set_model_parameter("lambda", lossRate);
  model->set_model_parameter("delta_avg", dupRate);
  model->set_model_parameter("sigma_hat", 1);
  model->calculate_EGb();
  auto res = model->p(ale);
  cout << "evalDated " << res << endl;
  return res;
}

double ALEEvaluation::evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
#ifdef UNDATED
  if (firstCall) {
    undatedDLModel.setInitialGeneTree(treeinfo);
  }
  return (double)std::log(undatedDLModel.pun(treeinfo));
#else
  char *newick = pll_utree_export_newick(treeinfo->root, 0);
  string geneTreeStr(newick);
  return evalDated(geneTreeStr,
    speciesTreeStr_,
    dupRate_,
    lossRate_,
    transferRate_);

#endif
}

