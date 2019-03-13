
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/DatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <cmath>

ReconciliationEvaluation::ReconciliationEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& map,
  RecModel recModel,
  bool rootedGeneTree):
  _model(recModel)
{
  reconciliationModel = getRecModelObject(recModel);
  reconciliationModel->init(speciesTree, map, rootedGeneTree);
}

void ReconciliationEvaluation::setRates(double dupRate, double lossRate,
  double transferRate)
{
  reconciliationModel->setRates(dupRate, lossRate, transferRate);
}

double ReconciliationEvaluation::evaluate(pll_utree_t *utree)
{
  return reconciliationModel->computeLogLikelihood(utree);
}

double ReconciliationEvaluation::evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  return evaluate(treeinfo->tree);
}

void ReconciliationEvaluation::invalidateCLV(int nodeIndex)
{
  reconciliationModel->invalidateCLV(nodeIndex);
}

shared_ptr<AbstractReconciliationModel> ReconciliationEvaluation::getRecModelObject(RecModel recModel)
{
  switch(recModel) {
  case UndatedDL:
    return  make_shared<UndatedDLModel>();
  case UndatedDTL:
    return  make_shared<UndatedDTLModel>();
  case DatedDL:
    return  make_shared<DatedDLModel>();
  }
  assert(false);
}

