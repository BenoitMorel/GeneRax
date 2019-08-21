
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <cmath>

ReconciliationEvaluation::ReconciliationEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& geneSpeciesMapping,
  RecModel recModel,
  bool rootedGeneTree):
  _model(recModel),
  _speciesCount(speciesTree->inner_count + speciesTree->tip_count)
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
  if (ratesVector.size() == 1) {
    dupRates = std::vector<double>(_speciesCount, ratesVector.getRates(0).rates[0]);
    lossRates = std::vector<double>(_speciesCount, ratesVector.getRates(0).rates[1]);
    transferRates = std::vector<double>(_speciesCount, ratesVector.getRates(0).rates[2]);
  } else {
    assert(ratesVector.size() == _speciesCount);
    for (auto &r: ratesVector.getRatesVector()) {
      dupRates.push_back(r.rates[0]); 
      lossRates.push_back(r.rates[1]); 
      transferRates.push_back(r.rates[2]); 
    }
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

std::shared_ptr<AbstractReconciliationModel<ScaledValue> > ReconciliationEvaluation::getRecModelObject(RecModel recModel)
{
  switch(recModel) {
  case UndatedDL:
    return  std::make_shared<UndatedDLModel<ScaledValue> >();
  case UndatedDTL:
    return  std::make_shared<UndatedDTLModel<ScaledValue> >();
  }
  assert(false);
  return 0;
}

