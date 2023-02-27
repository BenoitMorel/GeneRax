#include "AleOptimizer.hpp" 
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <util/Paths.hpp>
#include <maths/Random.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <search/DatedSpeciesTreeSearch.hpp>


static std::shared_ptr<MultiModel> createModel(SpeciesTree &speciesTree,
  const FamilyInfo &family,
  const ModelParameters &modelParameters,
  const std::vector<Highway> &highways,
  bool highPrecision)
{
  std::shared_ptr<MultiModel> model;
  GeneSpeciesMapping mapping;
  mapping.fill(family.mappingFile, family.startingGeneTree);
  const auto &info = modelParameters.info;  
  switch (info.model) {
  case RecModel::UndatedDL:
    if (highPrecision) {
      model = std::make_shared<UndatedDLMultiModel<ScaledValue> >(
        speciesTree.getTree(),
        mapping,
        info,
        family.ccp);
    } else {
      model = std::make_shared<UndatedDLMultiModel<double> >(
        speciesTree.getTree(),
        mapping,
        info,
        family.ccp);
    }
    break;
  case RecModel::UndatedDTL:
    if (highPrecision) {
      model = std::make_shared<UndatedDTLMultiModel<ScaledValue> >(
        speciesTree.getDatedTree(),
        mapping,
        info,
        family.ccp);
    } else {
      model = std::make_shared<UndatedDTLMultiModel<double> >(
        speciesTree.getDatedTree(),
        mapping,
        info,
        family.ccp);
    }
    break;
  default:
    assert(false);
    break;
  }
  std::vector<std::vector<double> > rates;
  unsigned int N = speciesTree.getTree().getNodesNumber();
  for (unsigned int d = 0; d < modelParameters.rates.dimensions(); ++d) {
    rates.push_back(std::vector<double>(N, modelParameters.rates[d]));
  }
  model->setRates(rates);
  model->setHighways(highways);
  return model;
}

GTSpeciesTreeLikelihoodEvaluator::GTSpeciesTreeLikelihoodEvaluator(
    SpeciesTree &speciesTree,
    ModelParameters &modelRates, 
    const Families &families,
    PerCoreGeneTrees &geneTrees):
  _speciesTree(speciesTree),
  _modelRates(modelRates),
  _families(families),
  _geneTrees(geneTrees)
  
{
  Logger::timed << "Initializing ccps and evaluators..." << std::endl;
  for (const auto &geneTree: _geneTrees.getTrees()) {
    auto &family = families[geneTree.familyIndex];
    _evaluations.push_back(createModel(_speciesTree,
          family,
          _modelRates,
          _highways,
          false));
    _highPrecisions.push_back(-1);
  }
  ParallelContext::barrier();
  unsigned int cladesNumber = 0;
  unsigned int worstFamily = 0;
  for (auto &evaluation: _evaluations) {
    cladesNumber += evaluation->getCCP().getCladesNumber();
    worstFamily = std::max(worstFamily, evaluation->getCCP().getCladesNumber());
  }
  unsigned int totalCladesNumber = cladesNumber;
  unsigned int maxCladesNumber = cladesNumber;
  ParallelContext::maxUInt(worstFamily);
  ParallelContext::maxUInt(maxCladesNumber);
  ParallelContext::sumUInt(totalCladesNumber);
  double averageCladesNumber = double(totalCladesNumber) / double(ParallelContext::getSize());
  Logger::timed << "Initializing ccps finished" << std::endl;
  Logger::timed << "Total number of clades: " << totalCladesNumber << std::endl;
  Logger::timed << "Load balancing: " << averageCladesNumber / maxCladesNumber << std::endl;
}

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihood()
{
  return computeLikelihoodFast();
}

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihoodFast()
{
  double sumLL = 0.0;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    auto ll = _evaluations[i]->computeLogLikelihood();
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto &family = _families[famIndex];
    if (_highPrecisions[i] == -1 && !std::isnormal(ll)) {
      // we are in low precision mode (we use double)
      // and it's not accurate enough, switch to
      // high precision mode
      _evaluations[i] = createModel(_speciesTree, 
          family,
          _modelRates,
          _highways,
          true);
      _highPrecisions[i] = 0;
      ll = _evaluations[i]->computeLogLikelihood();
    }
    if (!std::isnormal(ll)) {
      std::cerr << "Error: ll=" << ll << " for family " << family.name << std::endl;
    }
    assert(std::isnormal(ll));
    if (_highPrecisions[i] >= 0 && _highPrecisions[i] % 20 == 0) {
      // we are in high precision mode, we now check if we can
      // switch to low precision mode to make computations faster
      auto ev = createModel(_speciesTree, 
          family,
          _modelRates,
          _highways,
          false);
      auto newLL = ev->computeLogLikelihood();
      if (std::isnormal(newLL)) {
        _evaluations[i] = ev;
        _highPrecisions[i] = -1;
      }
    }
    if (_highPrecisions[i] >= 0) { 
      _highPrecisions[i]++;
    }
    sumLL += ll;
  }
  printHightPrecisionCount();
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

void GTSpeciesTreeLikelihoodEvaluator::setAlpha(double alpha)
{
  for (auto evaluation: _evaluations) {
    evaluation->setAlpha(alpha);
  }
}
  
void GTSpeciesTreeLikelihoodEvaluator::printHightPrecisionCount()
{
  unsigned int high = 0;
  unsigned int low = 0;
  for (auto v : _highPrecisions) {
    if (v >= 0) {
      high++;
    } else {
      low++;
    }
  }
  ParallelContext::sumUInt(high);
  ParallelContext::sumUInt(low);
  //Logger::info << " Double: " << low << " scaledvalue: " << high << std::endl;
}

void GTSpeciesTreeLikelihoodEvaluator::onSpeciesTreeChange(
    const std::unordered_set<corax_rnode_t *> *nodesToInvalidate)
{
  for (auto &evaluation: _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}

void GTSpeciesTreeLikelihoodEvaluator::onSpeciesDatesChange()
{
  for (auto &evaluation: _evaluations) {
    evaluation->onSpeciesDatesChange();
  }
}

void GTSpeciesTreeLikelihoodEvaluator::addHighway(const Highway &highway)
{
  _highways.push_back(highway);
  for (auto &evaluation: _evaluations) {
    evaluation->setHighways(_highways);
  }
}

void GTSpeciesTreeLikelihoodEvaluator::removeHighway()
{
  _highways.pop_back();
  for (auto &evaluation: _evaluations) {
    evaluation->setHighways(_highways);
  }
}

class GTEvaluatorFunction: public FunctionToOptimize
{
public:
  GTEvaluatorFunction(ReconciliationModelInterface &evaluation,
      RecModelInfo &info,
      unsigned int speciesNodeNumber):
    _evaluation(evaluation),
    _info(info),
    _speciesNodeNumber(speciesNodeNumber){}

  virtual double evaluate(Parameters &parameters) {
    parameters.ensurePositivity();
    unsigned int freeParameters = Enums::freeParameters(_info.model);
    if (!freeParameters) {
      return _evaluation.computeLogLikelihood();
    }
    assert(parameters.dimensions());
    assert(0 == parameters.dimensions() % freeParameters);
    std::vector<std::vector<double> > rates;
    rates.resize(freeParameters);
    for (auto &r: rates) {
      r.resize(_speciesNodeNumber);
    }
    // this handles both per-species and global rates
    for (unsigned int d = 0; d < rates.size(); ++d) {
      for (unsigned int e = 0; e < _speciesNodeNumber; ++e) {
        (rates[d])[e] = parameters[(e * rates.size() + d) % parameters.dimensions()];
      }
    }
    _evaluation.setRates(rates);
    
    double res = _evaluation.computeLogLikelihood();
    parameters.setScore(res);
    return res;
  }

private:
  ReconciliationModelInterface &_evaluation;
  RecModelInfo &_info;
  unsigned int _speciesNodeNumber;
};


class GTMultiEvaluatorsFunction: public FunctionToOptimize
{
public:
  GTMultiEvaluatorsFunction(PerCoreMultiEvaluation &evaluations,
      RecModelInfo &info, 
      unsigned int speciesNodeNumber)
  {
    for (auto &evaluation: evaluations) {
      _functions.push_back(std::make_shared<GTEvaluatorFunction>(*evaluation, 
            info,
            speciesNodeNumber));
    }
  }

  virtual double evaluate(Parameters &parameters) {
    double res = 0.0;
    for (auto &function: _functions) {
      res += function->evaluate(parameters);
    }
    ParallelContext::sumDouble(res);
    parameters.setScore(res);
    return res;
  }
private:
  std::vector<std::shared_ptr<GTEvaluatorFunction>> _functions;
};

double GTSpeciesTreeLikelihoodEvaluator::optimizeModelRates(bool thorough)
{
  OptimizationSettings settings;
  double ll = computeLikelihood();
  Logger::timed << "[SpeciesSearch] Optimizing model rates ";
  if (thorough) {
    Logger::info << "(thorough)" << std::endl;
  } else {
    Logger::info << "(light)" << std::endl;
  }
  if (!thorough) {
    settings.lineSearchMinImprovement = 10.0;
    settings.minAlpha = 0.01;
    settings.optimizationMinImprovement = std::max(3.0, ll / 1000.0);
  }
  GTMultiEvaluatorsFunction function(_evaluations, 
      _modelRates.info, 
      _speciesTree.getTree().getNodesNumber());
  _modelRates.rates = DTLOptimizer::optimizeParameters(
      function, 
      _modelRates.rates, 
      settings);
  ll = computeLikelihood();
  Logger::timed << "[Species search]   After model rate opt, ll=" << ll << " rates: " << _modelRates.rates << std::endl;
  ll = optimizeGammaRates();
  ll = computeLikelihood();
  return ll;
}

static double callback(void *p, double x)
{
  auto *evaluator = (GTSpeciesTreeLikelihoodEvaluator *)p;
  evaluator->setAlpha(x);
  auto ll = evaluator->computeLikelihood();
  return -ll;
}

double GTSpeciesTreeLikelihoodEvaluator::optimizeGammaRates()
{
  auto gammaCategories = _modelRates.info.gammaCategories;
  auto ll = computeLikelihood();
  if (gammaCategories == 1) {
    return ll;
  }
  double minAlpha = CORAX_OPT_MIN_ALPHA;  
  double maxAlpha = CORAX_OPT_MAX_ALPHA;  
  double startingAlpha = 1.0;
  double tolerance = 0.1;
  double f2x = 1.0;
  double alpha = corax_opt_minimize_brent(minAlpha,
                                  startingAlpha,
                                  maxAlpha,
                                  tolerance,
                                  &ll,
                                  &f2x,
                                  (void *)this,
                                  &callback);
  setAlpha(alpha);
  std::vector<double> categories(_modelRates.info.gammaCategories);
  corax_compute_gamma_cats(alpha, categories.size(), &categories[0], 
      CORAX_GAMMA_RATES_MEAN);
  Logger::timed << "[Species search]   After gamma cat  opt, ll=" << ll << std::endl;
  Logger::info << "alpha = " << alpha << std::endl;
  Logger::info << "rate categories: ";
  for (auto c: categories) {
    Logger::info << c << " ";
  }
  Logger::info << std::endl;
  return ll;
}
  
void GTSpeciesTreeLikelihoodEvaluator::getTransferInformation(SpeciesTree &speciesTree,
    TransferFrequencies &transferFrequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers)
{
  // this is duplicated code from Routines...
  const auto labelToId = speciesTree.getTree().getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getTree().getDeterministicIdToLabel();
  const unsigned int labelsNumber = idToLabel.size();
  const VectorUint zeros(labelsNumber, 0);
  transferFrequencies.count = MatrixUint(labelsNumber, zeros);
  transferFrequencies.idToLabel = idToLabel;
  perSpeciesEvents = PerSpeciesEvents(speciesTree.getTree().getNodesNumber());
  auto infoCopy = _modelRates.info;
  infoCopy.originationStrategy = OriginationStrategy::UNIFORM;
  for (const auto &geneTree: _geneTrees.getTrees()) {
    auto &family = (_families)[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    UndatedDTLMultiModel<ScaledValue> evaluation(
        speciesTree.getDatedTree(),
        mapping,
        infoCopy,
        family.ccp);
    
    evaluation.computeLogLikelihood();
    // warning, this might make the random state 
    // inconsistent between the MPI ranks
    // ParallelContext::makeRandConsistent() needs to be called 
    // right after the loop
    Scenario scenario;
    evaluation.inferMLScenario(scenario, true);
    scenario.countTransfers(labelToId, 
        transferFrequencies.count);
    scenario.gatherReconciliationStatistics(perSpeciesEvents);
    potentialTransfers.addScenario(scenario);
  }
  ParallelContext::makeRandConsistent();
  for (unsigned int i = 0; i < labelsNumber; ++i) {
    ParallelContext::sumVectorUInt(transferFrequencies.count[i]); 
  }
  perSpeciesEvents.parallelSum();
  assert(ParallelContext::isRandConsistent());
}

void GTSpeciesTreeLikelihoodEvaluator::fillPerFamilyLikelihoods(
    PerFamLL &perFamLL)
{
  ParallelContext::barrier();
  std::vector<double> localLL;
  for (auto &evaluation: _evaluations) {
    localLL.push_back(evaluation->computeLogLikelihood());
  }
  ParallelContext::concatenateHetherogeneousDoubleVectors(
      localLL, perFamLL);
}


