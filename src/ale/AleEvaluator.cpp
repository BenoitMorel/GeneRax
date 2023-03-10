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
//  highPrecision = true;
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
    bool optimizeRates,
    const Families &families,
    PerCoreGeneTrees &geneTrees,
    const std::string &outputDir):
  _speciesTree(speciesTree),
  _modelRates(modelRates),
  _optimizeRates(optimizeRates),
  _families(families),
  _geneTrees(geneTrees),
  _highPrecisions(_geneTrees.getTrees().size(), -1),
  _outputDir(outputDir)
{
  Logger::timed << "Initializing ccps and evaluators..." << std::endl;
  _evaluations.resize(_geneTrees.getTrees().size());
  for (unsigned int i = 0; i < _geneTrees.getTrees().size(); ++i) {
    resetEvaluation(i, false);
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


void GTSpeciesTreeLikelihoodEvaluator::resetEvaluation(unsigned int i, bool highPrecision)
{
  auto famIndex = _geneTrees.getTrees()[i].familyIndex;
  auto &family = _families[famIndex];
  _evaluations[i] = createModel(_speciesTree, 
      family,
      _modelRates,
      _highways,
      highPrecision);
  _highPrecisions[i] = highPrecision;
  auto ll = _evaluations[i]->computeLogLikelihood();
  if (highPrecision) {
    _highPrecisions[i] = 1;
  } else {
    _highPrecisions[i] = -1;
    if (!std::isnormal(ll)) {
      resetEvaluation(i, true);
    }
  }
}

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihoodFast()
{
  double sumLL = 0.0;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = _evaluations[i]->computeLogLikelihood();
    auto &family = _families[famIndex];
    if (_highPrecisions[i] == -1 && !std::isnormal(ll)) {
      // we are in low precision mode (we use double)
      // and it's not accurate enough, switch to
      // high precision mode
      resetEvaluation(i, true);
      ll = _evaluations[i]->computeLogLikelihood();
    }
    if (!std::isnormal(ll)) {
      std::cerr << "Error: ll=" << ll << " for family " << family.name << std::endl;
    }
    assert(std::isnormal(ll));
    if (_highPrecisions[i] >= 0 && _highPrecisions[i] % 20 == 0) {
      // we are in high precision mode, we now check if we can
      // switch to low precision mode to make computations faster
        resetEvaluation(i, false);
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


class DTLParametersOptimizer: public FunctionToOptimize
{
public:
  DTLParametersOptimizer(GTSpeciesTreeLikelihoodEvaluator &evaluator):
    _evaluator(evaluator)
  {}

  virtual double evaluate(Parameters &parameters) {
    parameters.ensurePositivity();
    _evaluator.setParameters(parameters);
    auto res = _evaluator.computeLikelihood();
    parameters.setScore(res);
    return res;

  }
private:
  GTSpeciesTreeLikelihoodEvaluator &_evaluator;
};
  

void GTSpeciesTreeLikelihoodEvaluator::setParameters(Parameters &parameters)
{
  unsigned int freeParameters = Enums::freeParameters(_modelRates.info.model);
  if (!freeParameters) {
    return;
  }
  assert(parameters.dimensions());
  assert(0 == parameters.dimensions() % freeParameters);
  std::vector<std::vector<double> > rates;
  rates.resize(freeParameters);
  auto speciesNodeNumber = _speciesTree.getTree().getNodesNumber();
  for (auto &r: rates) {
    r.resize(speciesNodeNumber);
  }
  // this handles both per-species and global rates
  for (unsigned int d = 0; d < rates.size(); ++d) {
    for (unsigned int e = 0; e < speciesNodeNumber; ++e) {
      (rates[d])[e] = parameters[(e * rates.size() + d) % parameters.dimensions()];
    }
  }
  for (auto evaluation: _evaluations) { 
    evaluation->setRates(rates);
  }
   
}

double GTSpeciesTreeLikelihoodEvaluator::optimizeModelRates(bool thorough)
{
  double ll = 0.0;
  if (_optimizeRates) {
    OptimizationSettings settings;
    ll = computeLikelihood();
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
    DTLParametersOptimizer function(*this);
    _modelRates.rates = DTLOptimizer::optimizeParameters(
        function, 
        _modelRates.rates, 
        settings);
    ll = computeLikelihood();
    Logger::timed << "[Species search]   After model rate opt, ll=" << ll << " rates: " << _modelRates.rates << std::endl;
  }
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
    bool ok = evaluation.inferMLScenario(scenario, true);
    assert(ok);
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

void GTSpeciesTreeLikelihoodEvaluator::sampleScenarios(unsigned int family, unsigned int samples,
      std::vector<Scenario> &scenarios)
{
  assert(family < _evaluations.size());
  scenarios = std::vector<Scenario>(samples);
  getEvaluation(family).computeLogLikelihood();
  for (unsigned int s = 0; s < samples; ++s) {
    bool ok  = getEvaluation(family).inferMLScenario(scenarios[s], true);
    if (!ok) {
      assert(_highPrecisions[family] == -1);
      resetEvaluation(family, true);
      sampleScenarios(family, samples, scenarios);
      return;
    }
  }
}
  
void GTSpeciesTreeLikelihoodEvaluator::savePerFamilyLikelihoodDiff(const std::string &output) 
{
  std::vector<unsigned int> indices;
  std::vector<double> likelihoods;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    auto famIndex = _geneTrees.getTrees()[i].familyIndex;
    auto ll = _evaluations[i]->computeLogLikelihood();
    indices.push_back(famIndex);
    likelihoods.push_back(ll);
  }
  std::vector<unsigned int> allIndices;
  std::vector<double> allLikelihoods;
  ParallelContext::concatenateHetherogeneousDoubleVectors(likelihoods,
  allLikelihoods);
  ParallelContext::concatenateHetherogeneousUIntVectors(indices, allIndices);
  assert(allLikelihoods.size() == _snapshotPerFamilyLL.size());
  ParallelOfstream os(output);
  for (unsigned int i = 0; i < allLikelihoods.size(); ++i) {
    auto &family = _families[allIndices[i]];
    auto ll = allLikelihoods[i];
    os << family.name << " " << ll - _snapshotPerFamilyLL[i]<< std::endl;
  }
}

void GTSpeciesTreeLikelihoodEvaluator::saveSnapshotPerFamilyLL()
{
  std::vector<double> likelihoods;
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    auto ll = _evaluations[i]->computeLogLikelihood();
    likelihoods.push_back(ll);
  }
  ParallelContext::concatenateHetherogeneousDoubleVectors(likelihoods,
  _snapshotPerFamilyLL);
}


