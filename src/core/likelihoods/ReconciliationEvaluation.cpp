
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <cmath>

ReconciliationEvaluation::ReconciliationEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& geneSpeciesMapping,
  RecModel recModel,
  bool rootedGeneTree):
  _model(recModel)
{
  reconciliationModel = getRecModelObject(recModel);
  reconciliationModel->init(speciesTree, geneSpeciesMapping, rootedGeneTree);
}

void ReconciliationEvaluation::setRates(double dupRate, double lossRate,
  double transferRate)
{
  reconciliationModel->setRates(dupRate, lossRate, transferRate);
}
  
void ReconciliationEvaluation::setRates(const DTLRatesVector &ratesVector)
{
  std::vector<double> dupRates;
  std::vector<double> lossRates;
  std::vector<double> transferRates;
  for (auto &r: ratesVector.getRatesVector()) {
    dupRates.push_back(r.rates[0]); 
    lossRates.push_back(r.rates[1]); 
    transferRates.push_back(r.rates[2]); 
  }
  reconciliationModel->setRates(dupRates, lossRates, transferRates);
}

double ReconciliationEvaluation::evaluate(pll_utree_t *utree)
{
  return reconciliationModel->computeLogLikelihood(utree);
}

double ReconciliationEvaluation::evaluate(std::shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  return evaluate(treeinfo->tree);
}

void ReconciliationEvaluation::invalidateCLV(unsigned int nodeIndex)
{
  reconciliationModel->invalidateCLV(nodeIndex);
}

std::shared_ptr<AbstractReconciliationModel> ReconciliationEvaluation::getRecModelObject(RecModel recModel)
{
  switch(recModel) {
  case UndatedDL:
    return  std::make_shared<UndatedDLModel>();
  case UndatedDTL:
    return  std::make_shared<UndatedDTLModel>();
  }
  assert(false);
  return 0;
}

