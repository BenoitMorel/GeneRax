
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedIDTLModel.hpp>
#include <cmath>
#include <IO/FileSystem.hpp>
#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>

double log(ScaledValue v) 
{
  return v.getLogValue();
}


ReconciliationEvaluation::ReconciliationEvaluation(PLLRootedTree  &speciesTree,
  PLLUnrootedTree &initialGeneTree,
  const GeneSpeciesMapping& geneSpeciesMapping,
  RecModel recModel,
  bool rootedGeneTree, 
  bool pruneSpeciesTree):
    _speciesTree(speciesTree),
    _initialGeneTree(initialGeneTree),
    _geneSpeciesMapping(geneSpeciesMapping),
    _rootedGeneTree(rootedGeneTree),
    _pruneSpeciesTree(pruneSpeciesTree),
    _model(recModel),
    _infinitePrecision(true)
{
  _evaluators = buildRecModelObject(_model, _infinitePrecision);
}
  
ReconciliationEvaluation::~ReconciliationEvaluation()
{
  delete _evaluators;
}


void ReconciliationEvaluation::setRates(const Parameters &parameters)
{
  assert(parameters.dimensions());
  unsigned int freeParameters = Enums::freeParameters(_model);
  assert(0 == parameters.dimensions() % freeParameters);
  _rates.resize(freeParameters);
  for (auto &r: _rates) {
    r.resize(_speciesTree.getNodesNumber());
  }
  // this handles both per-species and global rates
  for (unsigned int d = 0; d < _rates.size(); ++d) {
    for (unsigned int e = 0; e < _speciesTree.getNodesNumber(); ++e) {
      (_rates[d])[e] = parameters[(e * _rates.size() + d) % parameters.dimensions()];
    }
  }
  _evaluators->setRates(_rates);
}
  
pll_unode_t *ReconciliationEvaluation::getRoot() 
{
  return _evaluators->getRoot();
}

void ReconciliationEvaluation::setRoot(pll_unode_t * root) 
{
  _evaluators->setRoot(root);
}

double ReconciliationEvaluation::evaluate(bool fastMode)
{
  fastMode = false;
  double res = _evaluators->computeLogLikelihood(fastMode);
  if (!_infinitePrecision && !std::isnormal(res)) {
    updatePrecision(true);  
    res = _evaluators->computeLogLikelihood(fastMode);
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
  _evaluators->invalidateCLV(nodeIndex);
}

void ReconciliationEvaluation::invalidateAllCLVs()
{
  _evaluators->invalidateAllCLVs();
}

void ReconciliationEvaluation::invalidateAllSpeciesCLVs()
{
  _evaluators->invalidateAllSpeciesCLVs();
}

ReconciliationModelInterface *ReconciliationEvaluation::buildRecModelObject(RecModel recModel, 
    bool infinitePrecision)
{
  ReconciliationModelInterface *res(nullptr);
  switch(recModel) {
  case RecModel::UndatedDL:
    if (infinitePrecision) {
      res = new UndatedDLModel<ScaledValue>(_speciesTree, _geneSpeciesMapping, _rootedGeneTree, _pruneSpeciesTree);
    } else {
      res = new UndatedDLModel<double>(_speciesTree, _geneSpeciesMapping, _rootedGeneTree, _pruneSpeciesTree);
    }
    break;
  case RecModel::UndatedDTL:
    if (infinitePrecision) {
      res = new UndatedDTLModel<ScaledValue>(_speciesTree, _geneSpeciesMapping, _rootedGeneTree, _pruneSpeciesTree);
    } else {
      res = new UndatedDTLModel<double>(_speciesTree, _geneSpeciesMapping, _rootedGeneTree, _pruneSpeciesTree);
    }
    break;
  case RecModel::UndatedIDTL:
    if (infinitePrecision) {
      res = new UndatedIDTLModel<ScaledValue>(_speciesTree, _geneSpeciesMapping, _rootedGeneTree, _pruneSpeciesTree);
    } else {
      res = new UndatedIDTLModel<double>(_speciesTree, _geneSpeciesMapping, _rootedGeneTree, _pruneSpeciesTree);
    }
    break;
  }
  res->setInitialGeneTree(_initialGeneTree.getRawPtr());
  return res;
}
  
void ReconciliationEvaluation::updatePrecision(bool infinitePrecision)
{
  if (infinitePrecision != _infinitePrecision) {
    _infinitePrecision = infinitePrecision;
    delete _evaluators;
    _evaluators = buildRecModelObject(_model, _infinitePrecision);
    _evaluators->setRates(_rates);
 }
}

void ReconciliationEvaluation::inferMLScenario(Scenario &scenario, bool stochastic) {
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate();
  assert(std::isfinite(ll) && ll < 0.0);
  _evaluators->inferMLScenario(scenario, stochastic);
  updatePrecision(infinitePrecision);
}
  
pll_unode_t *ReconciliationEvaluation::computeMLRoot() 
{
  return  _evaluators->computeMLRoot();
}
  
pll_unode_t *ReconciliationEvaluation::inferMLRoot()
{
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
  assert(_evaluators);
  _evaluators->onSpeciesTreeChange(nodesToInvalidate);
}

void ReconciliationEvaluation::setPartialLikelihoodMode(PartialLikelihoodMode mode) 
{ 
  _evaluators->setPartialLikelihoodMode(mode);
}
  
void ReconciliationEvaluation::rollbackToLastState() 
{
  _evaluators->rollbackToLastState();
}

