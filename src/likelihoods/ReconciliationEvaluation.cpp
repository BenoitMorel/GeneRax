
#include "ReconciliationEvaluation.hpp"
#include <Logger.hpp>
#include <likelihoods/reconciliation_models/DatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <cmath>

ReconciliationEvaluation::ReconciliationEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& map,
  Arguments::ReconciliationModel model)
{
  if (model == Arguments::UndatedDL) {
    reconciliationModel = make_shared<UndatedDLModel>(speciesTree, map);
  } else if (model == Arguments::DatedDL) {
    reconciliationModel = make_shared<DatedDLModel>(speciesTree, map);
  } else {
    Logger::error << "Invalid reconciliation model!" << endl;
    exit(1);
  }
}

void ReconciliationEvaluation::setRates(double dupRate, double lossRate,
  double transferRate)
{
  reconciliationModel->setRates(dupRate, lossRate, transferRate);
}

double ReconciliationEvaluation::evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  return reconciliationModel->computeLogLikelihood(treeinfo);
}

void ReconciliationEvaluation::invalidateCLV(int nodeIndex)
{
  reconciliationModel->invalidateCLV(nodeIndex);
}
