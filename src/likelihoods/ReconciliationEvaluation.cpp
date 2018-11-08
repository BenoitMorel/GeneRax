
#include "ReconciliationEvaluation.hpp"
#include <likelihoods/ale/UndatedDLModel.hpp>
#include <likelihoods/ale/UndatedDTLModel.hpp>
#include <cmath>

ReconciliationEvaluation::ReconciliationEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& map,
  bool transfers):
  firstCall(true)
{
  if (transfers) {
    reconciliationModel = make_shared<UndatedDTLModel>();
  } else {
    reconciliationModel = make_shared<UndatedDLModel>();
  }
}

void ReconciliationEvaluation::setRates(double dupRate, double lossRate,
  double transferRate)
{
  reconciliationModel->setRates(dupRate, lossRate);
}

double ReconciliationEvaluation::evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  if (firstCall) {
    reconciliationModel->setInitialGeneTree(treeinfo);
  }
  return (double)log(reconciliationModel->computeLikelihood(treeinfo));
}

