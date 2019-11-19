
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <cmath>
#include <IO/FileSystem.hpp>

double log(ScaledValue v) 
{
  return v.getLogValue();
}


ReconciliationEvaluation::ReconciliationEvaluation(PLLRootedTree  &speciesTree,
  PLLUnrootedTree &initialGeneTree,
  const GeneSpeciesMapping& geneSpeciesMapping,
  RecModel recModel,
  bool rootedGeneTree):
    _speciesTree(speciesTree),
    _initialGeneTree(initialGeneTree),
    _geneSpeciesMapping(geneSpeciesMapping),
    _rootedGeneTree(rootedGeneTree),
    _model(recModel),
    _infinitePrecision(true)
{
  _reconciliationModel = buildRecModelObject(_model, _infinitePrecision);
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
  _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates);
}

double ReconciliationEvaluation::evaluate(bool fastMode)
{
  double res = _reconciliationModel->computeLogLikelihood(fastMode);
  if (!_infinitePrecision && !std::isnormal(res)) {
    Logger::info << "not normal ll: " << res << std::endl;
    updatePrecision(true);  
    res = _reconciliationModel->computeLogLikelihood(fastMode);
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

void ReconciliationEvaluation::invalidateAllCLVs()
{
  _reconciliationModel->invalidateAllCLVs();
}

void ReconciliationEvaluation::invalidateAllSpeciesCLVs()
{
  _reconciliationModel->invalidateAllSpeciesCLVs();
}

std::unique_ptr<ReconciliationModelInterface> ReconciliationEvaluation::buildRecModelObject(RecModel recModel, 
    bool infinitePrecision)
{
  std::unique_ptr<ReconciliationModelInterface> res(nullptr);
  switch(recModel) {
  case RecModel::UndatedDL:
    if (infinitePrecision) {
      res = std::make_unique<UndatedDLModel<ScaledValue> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    } else {
      res = std::make_unique<UndatedDLModel<double> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    }
    break;
  case RecModel::UndatedDTL:
    if (infinitePrecision) {
      res = std::make_unique<UndatedDTLModel<ScaledValue> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    } else {
      res = std::make_unique<UndatedDTLModel<double> >(_speciesTree, _geneSpeciesMapping, _rootedGeneTree);
    }
    break;
  }
  res->setInitialGeneTree(_initialGeneTree.getRawPtr());
  return res;
}
  
void ReconciliationEvaluation::updatePrecision(bool infinitePrecision)
{
  if (infinitePrecision != _infinitePrecision) {
    Logger::info << "UPDATE PRECISION " << infinitePrecision << std::endl;
    _infinitePrecision = infinitePrecision;
    _reconciliationModel = buildRecModelObject(_model, _infinitePrecision);
    _reconciliationModel->setRates(_dupRates, _lossRates, _transferRates);
 }
}

void ReconciliationEvaluation::inferMLScenario(Scenario &scenario) {
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate();
  assert(std::isfinite(ll) && ll < 0.0);
  _reconciliationModel->inferMLScenario(scenario);
  updatePrecision(infinitePrecision);
}
  
pll_unode_t *ReconciliationEvaluation::computeMLRoot() 
{
  return  _reconciliationModel->computeMLRoot();
}
  
pll_unode_t *ReconciliationEvaluation::inferMLRoot()
{
  Logger::info << "infer ML root" << std::endl;
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate(); 
  assert(std::isfinite(ll) && ll < 0.0);
  auto res = computeMLRoot();
  updatePrecision(infinitePrecision);
  assert(res);
  return res;
}

void ReconciliationEvaluation::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  assert(_reconciliationModel);
  _reconciliationModel->onSpeciesTreeChange(nodesToInvalidate);
}
