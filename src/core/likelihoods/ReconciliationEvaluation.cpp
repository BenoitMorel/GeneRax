
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


ReconciliationEvaluation::ReconciliationEvaluation(PLLRootedTree  &speciesTree,
  const GeneSpeciesMapping& geneSpeciesMapping,
  RecModel recModel,
  bool rootedGeneTree):
    _speciesTree(speciesTree),
    _geneSpeciesMapping(geneSpeciesMapping),
    _rootedGeneTree(rootedGeneTree),
    _model(recModel),
    _infinitePrecision(false)
{
  _reconciliationModel = buildRecModelObject(_model, _infinitePrecision);
}


void ReconciliationEvaluation::setTransferFrequencies(const Parameters &parameters)
{
  assert(parameters.dimensions() == _speciesTree.getNodesNumber() * _speciesTree.getNodesNumber());
  _transferFrequencies.resize(_speciesTree.getNodesNumber());
  unsigned int index = 0;
  for (unsigned int e = 0; e < _speciesTree.getNodesNumber(); ++e) {
    _transferFrequencies[e].resize(_speciesTree.getNodesNumber());
    for (unsigned int h = 0; h < _speciesTree.getNodesNumber(); ++h) {
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
    r->resize(_speciesTree.getNodesNumber());
  }
  // this handles both per-species and global rates
  for (unsigned int d = 0; d < rates.size(); ++d) {
    for (unsigned int e = 0; e < _speciesTree.getNodesNumber(); ++e) {
      (*rates[d])[e] = parameters[(e * rates.size() + d) % parameters.dimensions()];
    }
  }
  if (_model == RecModel::UndatedDTLAdvanced) {
    std::string transferFrequenciesFile = "/tmp/transferFrequencies.txt";
    if (FileSystem::exists(transferFrequenciesFile)) {
      Parameters hackParameters;
      hackParameters.load(transferFrequenciesFile);
      setTransferFrequencies(hackParameters);
    } else {
      _transferFrequencies.resize(_speciesTree.getNodesNumber());
      for (unsigned int e = 0; e < _speciesTree.getNodesNumber(); ++e) {
        _transferFrequencies[e] = std::vector<double>(_speciesTree.getNodesNumber(), 
            1.0 / _speciesTree.getNodesNumber());
      }
    }
  }
  _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates, _transferFrequencies);
}

double ReconciliationEvaluation::evaluate(PLLUnrootedTree &geneTree)
{
  auto utree = geneTree.getRawPtr();
  double res = _reconciliationModel->computeLogLikelihood(utree);
  if (!_infinitePrecision && !std::isnormal(res)) {
    updatePrecision(true);  
    res = _reconciliationModel->computeLogLikelihood(utree);
    updatePrecision(false);  
  }
  if (!std::isnormal(res)) {
    std::cerr << "wrong reconciliation ll " << res << std::endl;
  }
  assert(std::isnormal(res));
  return res;
}

void ReconciliationEvaluation::invalidateCLV(unsigned int nodeIndex)
{
  _reconciliationModel->invalidateCLV(nodeIndex);
}

std::unique_ptr<ReconciliationModelInterface> ReconciliationEvaluation::buildRecModelObject(RecModel recModel, 
    bool infinitePrecision)
{
  switch(recModel) {
  case RecModel::UndatedDL:
    if (infinitePrecision) {
      return  std::make_unique<UndatedDLModel<ScaledValue> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    } else {
      return  std::make_unique<UndatedDLModel<double> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    }
  case RecModel::UndatedDTL:
    if (infinitePrecision) {
      return  std::make_unique<UndatedDTLModel<ScaledValue> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    } else {
      return  std::make_unique<UndatedDTLModel<double> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    }
  case RecModel::UndatedDTLAdvanced:
    if (infinitePrecision) {
      return  std::make_unique<UndatedDTLModelAdvanced<ScaledValue> >(_speciesTree, 
          _geneSpeciesMapping, 
          _rootedGeneTree);
    } else {
      return  std::make_unique<UndatedDTLModelAdvanced<double> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
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
    _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates, _transferFrequencies);
 }
}

void ReconciliationEvaluation::inferMLScenario(PLLUnrootedTree &geneTree, Scenario &scenario) {
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate(geneTree);
  assert(std::isfinite(ll) && ll < 0.0);
  _reconciliationModel->inferMLScenario(scenario);
  updatePrecision(infinitePrecision);
}
  
pll_unode_t *ReconciliationEvaluation::computeMLRoot() 
{
  return  _reconciliationModel->computeMLRoot();
}
  
pll_unode_t *ReconciliationEvaluation::inferMLRoot(PLLUnrootedTree &geneTree)
{
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate(geneTree); 
  assert(std::isfinite(ll) && ll < 0.0);
  auto res = computeMLRoot();
  updatePrecision(infinitePrecision);
  assert(res);
  return res;
}

