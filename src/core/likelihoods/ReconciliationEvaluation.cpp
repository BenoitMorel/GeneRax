
#include "ReconciliationEvaluation.hpp"
#include <IO/Logger.hpp>
#include <likelihoods/reconciliation_models/SimpleDSModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDLModel.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <likelihoods/reconciliation_models/ParsimonyDModel.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>
#include <cmath>
#include <IO/FileSystem.hpp>

double log(ScaledValue v) 
{
  return v.getLogValue();
}


ReconciliationEvaluation::ReconciliationEvaluation(PLLRootedTree  &speciesTree,
  PLLUnrootedTree &initialGeneTree,
  const GeneSpeciesMapping& geneSpeciesMapping,
  const RecModelInfo &recModelInfo):
    _speciesTree(speciesTree),
    _initialGeneTree(initialGeneTree),
    _geneSpeciesMapping(geneSpeciesMapping),
    _recModelInfo(recModelInfo),
    _infinitePrecision(true)
{
  _evaluators = buildRecModelObject(_recModelInfo.model, 
      _infinitePrecision);
}
  
ReconciliationEvaluation::~ReconciliationEvaluation()
{
  delete _evaluators;
}


void ReconciliationEvaluation::setRates(const Parameters &parameters)
{
  unsigned int freeParameters = Enums::freeParameters(_recModelInfo.model);
  if (!freeParameters) {
    return;
  }
  assert(parameters.dimensions());
  assert(0 == parameters.dimensions() % freeParameters);
  _rates.resize(freeParameters);
  for (auto &r: _rates) {
    r.resize(_speciesTree.getNodeNumber());
  }
  // this handles both per-species and global rates
  for (unsigned int d = 0; d < _rates.size(); ++d) {
    for (unsigned int e = 0; e < _speciesTree.getNodeNumber(); ++e) {
      (_rates[d])[e] = parameters[(e * _rates.size() + d) % parameters.dimensions()];
    }
  }
  _evaluators->setRates(_rates);
}
  
corax_unode_t *ReconciliationEvaluation::getRoot() 
{
  return _evaluators->getRoot();
}

void ReconciliationEvaluation::setRoot(corax_unode_t * root) 
{
  _evaluators->setRoot(root);
}

double ReconciliationEvaluation::evaluate()
{
  return _evaluators->computeLogLikelihood();
}
  
void ReconciliationEvaluation::enableMADRooting(bool enable)
{
  _evaluators->enableMADRooting(enable);
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

GTBaseReconciliationInterface *ReconciliationEvaluation::buildRecModelObject(RecModel recModel, 
    bool infinitePrecision)
{
  GTBaseReconciliationInterface *res(nullptr);
  switch(recModel) {
  case RecModel::UndatedDL:
    if (infinitePrecision) {
      res = new UndatedDLModel<ScaledValue>(_speciesTree, _geneSpeciesMapping, _recModelInfo);
    } else {
      res = new UndatedDLModel<double>(_speciesTree, _geneSpeciesMapping, _recModelInfo);
    }
    break;
  case RecModel::UndatedDTL:
    if (infinitePrecision) {
      res = new UndatedDTLModel<ScaledValue>(_speciesTree, _geneSpeciesMapping, _recModelInfo);
    } else {
      res = new UndatedDTLModel<double>(_speciesTree, _geneSpeciesMapping, _recModelInfo);
    }
    break;
  case RecModel::ParsimonyD:
    res = new ParsimonyDModel(_speciesTree, _geneSpeciesMapping,
        _recModelInfo);
    break;
  case RecModel::SimpleDS:
    if (infinitePrecision) {
      res = new SimpleDSModel<ScaledValue>(_speciesTree, _geneSpeciesMapping, _recModelInfo);
    } else {
      res = new SimpleDSModel<double>(_speciesTree, _geneSpeciesMapping, _recModelInfo);
    }
    break;
  }
  res->setInitialGeneTree(_initialGeneTree);
  return res;
}
  
void ReconciliationEvaluation::updatePrecision(bool infinitePrecision)
{
  if (infinitePrecision != _infinitePrecision) {
    _infinitePrecision = infinitePrecision;
    delete _evaluators;
    _evaluators = buildRecModelObject(_recModelInfo.model, 
      _infinitePrecision);
    _evaluators->setRates(_rates);
 }
}
void ReconciliationEvaluation::inferMLScenario(Scenario &scenario, bool stochastic) {
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate();
  assert(std::isfinite(ll) && ll <= 0.0);
  _evaluators->inferMLScenario(scenario, stochastic);
  updatePrecision(infinitePrecision);
}
  
corax_unode_t *ReconciliationEvaluation::computeMLRoot() 
{
  return  _evaluators->computeMLRoot();
}
  
corax_unode_t *ReconciliationEvaluation::inferMLRoot()
{
  auto infinitePrecision = _infinitePrecision;
  updatePrecision(true);
  auto ll = evaluate(); 
  assert(std::isfinite(ll) && ll <= 0.0);
  auto res = computeMLRoot();
  updatePrecision(infinitePrecision);
  assert(res);
  return res;
}

void ReconciliationEvaluation::onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate)
{
  assert(_evaluators);
  _evaluators->onSpeciesTreeChange(nodesToInvalidate);
}

void ReconciliationEvaluation::setPartialLikelihoodMode(PartialLikelihoodMode mode) 
{ 
  _evaluators->setPartialLikelihoodMode(mode);
}
  
