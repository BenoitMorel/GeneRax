
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModelAdvanced.hpp>
#include <cmath>
#include <IO/FileSystem.hpp>

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


void ReconciliationEvaluation::setTransferFrequencies(const Parameters &parameters)
{
  assert(parameters.dimensions() == _speciesCount * _speciesCount);
  _transferFrequencies.resize(_speciesCount);
  unsigned int index = 0;
  for (unsigned int e = 0; e < _speciesCount; ++e) {
    _transferFrequencies[e].resize(_speciesCount);
    for (unsigned int h = 0; h < _speciesCount; ++h) {
      _transferFrequencies[e][h] = parameters[index++];
      assert(std::isfinite(_transferFrequencies[e][h]));
    }
  }
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
  if (_model == UndatedDTLAdvanced) {
    std::string transferFrequenciesFile = "/tmp/transferFrequencies.txt";
    if (FileSystem::exists(transferFrequenciesFile)) {
      Parameters parameters;
      parameters.load(transferFrequenciesFile);
      setTransferFrequencies(parameters);
    } else {
      _transferFrequencies.resize(_speciesCount);
      for (unsigned int e = 0; e < _speciesCount; ++e) {
        _transferFrequencies[e] = std::vector<double>(_speciesCount, 1.0 / _speciesCount);
      }
    }
  }
  _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates, _transferFrequencies);
}

double ReconciliationEvaluation::evaluate(pll_utree_t *utree)
{
  double res = _reconciliationModel->computeLogLikelihood(utree);
  if (!_infinitePrecision && !std::isfinite(res)) {
    updatePrecision(true);  
    res = _reconciliationModel->computeLogLikelihood(utree);
    updatePrecision(false);  
  }
  return res;
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
    _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates, _transferFrequencies);
 }
}
void ReconciliationEvaluation::inferMLScenario(pll_utree_t *tree, Scenario &scenario) {
  assert(tree);
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate(tree);
  assert(std::isfinite(ll) && ll < 0.0);
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
  assert(std::isfinite(ll) && ll < 0.0);
  auto res = computeMLRoot();
  updatePrecision(infinitePrecision);
  assert(res);
  return res;
}

