
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/DatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <cmath>

ReconciliationEvaluation::ReconciliationEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& map,
  Arguments::ReconciliationModel model)
{
  if (model == Arguments::UndatedDL) {
    reconciliationModel = make_shared<UndatedDLModel>();
  } else if (model == Arguments::UndatedDTL) {
    reconciliationModel = make_shared<UndatedDTLModel>();
  } else if (model == Arguments::DatedDL) {
    reconciliationModel = make_shared<DatedDLModel>();
  } else {
    Logger::error << "Invalid reconciliation model!" << endl;
    exit(1);
  }
  reconciliationModel->init(speciesTree, map);
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
