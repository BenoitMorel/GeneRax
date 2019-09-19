
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModelAdvanced.hpp>
#include <cmath>

double log(ScaledValue v) 
{
  return v.getLogValue();
}

ReconciliationEvaluation::ReconciliationEvaluation(pll_rtree_t *speciesTree,
  const GeneSpeciesMapping& geneSpeciesMapping,
  RecModel recModel,
  bool rootedGeneTree):
    _speciesTree(speciesTree),
    _speciesCount(speciesTree->inner_count + speciesTree->tip_count),
    _geneSpeciesMapping(geneSpeciesMapping),
    _rootedGeneTree(rootedGeneTree),
    _model(recModel),
    _infinitePrecision(false)
{
  _reconciliationModel = buildRecModelObject(_model, _infinitePrecision);
  _reconciliationModel->init(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
}


void ReconciliationEvaluation::setRates(const Parameters &parameters)
{
  assert(parameters.dimensions());
  std::vector<std::vector<double> *> rates;
  rates.push_back(&_dupRates);
  rates.push_back(&_lossRates);
  if (implementsTransfers()) {
    rates.push_back(&_transferRates);
  }
  for (auto r: rates) {
    r->resize(_speciesCount);
  }
  // this handles both per-species and global rates
  for (unsigned int d = 0; d < rates.size(); ++d) {
    for (unsigned int e = 0; e < _speciesCount; ++e) {
      (*rates[d])[e] = parameters[(e * rates.size() + d) % parameters.dimensions()];
    }
  }
  _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates);
}

void ReconciliationEvaluation::setRates(double dupRate, double lossRate,
  double transferRate)
{
  DTLRatesVector v(DTLRates(dupRate, lossRate, transferRate));;
  setRates(v);
}
  
void ReconciliationEvaluation::setRates(const DTLRatesVector &ratesVector)
{
  _dupRates.clear();
  _lossRates.clear();
  _transferRates.clear();
  if (ratesVector.size() == 1) {
    _dupRates = std::vector<double>(_speciesCount, ratesVector.getRates(0).rates[0]);
    _lossRates = std::vector<double>(_speciesCount, ratesVector.getRates(0).rates[1]);
    _transferRates = std::vector<double>(_speciesCount, ratesVector.getRates(0).rates[2]);
  } else {
    assert(ratesVector.size() == _speciesCount);
    for (auto &r: ratesVector.getRatesVector()) {
      _dupRates.push_back(r.rates[0]); 
      _lossRates.push_back(r.rates[1]); 
      _transferRates.push_back(r.rates[2]); 
    }
  }
  _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates);
}

double ReconciliationEvaluation::evaluate(pll_utree_t *utree)
{
  double res = _reconciliationModel->computeLogLikelihood(utree);
  if (!_infinitePrecision && !std::isnormal(res)) {
    updatePrecision(true);  
    res = _reconciliationModel->computeLogLikelihood(utree);
    updatePrecision(false);  
  }
  return res;
}

double ReconciliationEvaluation::evaluate(std::shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  return evaluate(treeinfo->tree);
}

void ReconciliationEvaluation::invalidateCLV(unsigned int nodeIndex)
{
  _reconciliationModel->invalidateCLV(nodeIndex);
}

std::unique_ptr<ReconciliationModelInterface> ReconciliationEvaluation::buildRecModelObject(RecModel recModel, bool infinitePrecision)
{
  switch(recModel) {
  case UndatedDL:
    if (infinitePrecision) {
      return  std::make_unique<UndatedDLModel<ScaledValue> >();
    } else {
      return  std::make_unique<UndatedDLModel<double> >();
    }
  case UndatedDTL:
    if (infinitePrecision) {
      return  std::make_unique<UndatedDTLModel<ScaledValue> >();
    } else {
      return  std::make_unique<UndatedDTLModel<double> >();
    }
  case UndatedDTLAdvanced:
    if (infinitePrecision) {
      return  std::make_unique<UndatedDTLModelAdvanced<ScaledValue> >();
    } else {
      return  std::make_unique<UndatedDTLModelAdvanced<double> >();
    }
  }
  assert(false);
  return 0;
}
  
void ReconciliationEvaluation::updatePrecision(bool infinitePrecision)
{
  if (infinitePrecision != _infinitePrecision) {
    _infinitePrecision = infinitePrecision;
    _reconciliationModel = buildRecModelObject(_model, _infinitePrecision);
    _reconciliationModel->init(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates);
  }
}
void ReconciliationEvaluation::inferMLScenario(pll_utree_t *tree, Scenario &scenario) {
  assert(tree);
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate(tree);
  assert(std::isnormal(ll) && ll < 0.0);
  _reconciliationModel->inferMLScenario(scenario);
  updatePrecision(infinitePrecision);
}
  
pll_unode_t *ReconciliationEvaluation::computeMLRoot() 
{
  return  _reconciliationModel->computeMLRoot();
}
  
pll_unode_t *ReconciliationEvaluation::inferMLRoot(pll_utree_t *tree)
{
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate(tree); 
  assert(std::isnormal(ll) && ll < 0.0);
  auto res = computeMLRoot();
  updatePrecision(infinitePrecision);
  assert(res);
  return res;
}

