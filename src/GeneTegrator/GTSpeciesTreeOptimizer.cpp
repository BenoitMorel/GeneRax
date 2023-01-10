#include "GTSpeciesTreeOptimizer.hpp" 
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <util/Paths.hpp>
#include <maths/Random.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <search/DatedSpeciesTreeSearch.hpp>

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

static std::shared_ptr<MultiModel> createModel(SpeciesTree &speciesTree,
  const FamilyInfo &family,
  const ModelParameters &modelParameters,
  bool highPrecision,
  int maxSamples = -1)
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
        family.startingGeneTree,
        maxSamples);
    } else {
      model = std::make_shared<UndatedDLMultiModel<double> >(
        speciesTree.getTree(),
        mapping,
        info,
        family.startingGeneTree,
        maxSamples);
    }
    break;
  case RecModel::UndatedDTL:
    if (highPrecision) {
      model = std::make_shared<UndatedDTLMultiModel<ScaledValue> >(
        speciesTree.getDatedTree(),
        mapping,
        info,
        family.startingGeneTree,
        maxSamples);
    } else {
      model = std::make_shared<UndatedDTLMultiModel<double> >(
        speciesTree.getDatedTree(),
        mapping,
        info,
        family.startingGeneTree,
        maxSamples);
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
          true);
      _highPrecisions[i] = 0;
      ll = _evaluations[i]->computeLogLikelihood();
      assert(std::isnormal(ll));
    }

    if (_highPrecisions[i] >= 0 && _highPrecisions[i] % 20 == 0) {
      // we are in high precision mode, we now check if we can
      // switch to low precision mode to make computations faster
      auto ev = createModel(_speciesTree, 
          family,
          _modelRates,
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
  Logger::timed << "[Species search]   After model rate opt, ll=" << ll << " initial rates: " << _modelRates.rates << std::endl;
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
    PerSpeciesEvents &perSpeciesEvents)
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
        family.startingGeneTree);
    
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




GTSpeciesTreeOptimizer::GTSpeciesTreeOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    const RecModelInfo &info,
    const std::string &outputDir):
  _speciesTree(std::make_unique<SpeciesTree>(speciesTreeFile)),
  _geneTrees(families, false, true),
  _info(info),
  _outputDir(outputDir),
  _searchState(*_speciesTree, 
      Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"))
{
  Parameters startingRates;
  switch (info.model) {
  case RecModel::UndatedDL:
    startingRates = Parameters(0.2, 0.2);
    break;
  case RecModel::UndatedDTL:
    startingRates = Parameters(0.2, 0.2, 0.1);
    break;
  default:
    assert(false);
    break;
  }
  _modelRates = ModelParameters(startingRates, 
      1, //_geneTrees.getTrees().size(), 
      info);
  _speciesTree->addListener(this);
  ParallelContext::barrier();
  _evaluator = std::make_unique<GTSpeciesTreeLikelihoodEvaluator>(
      *_speciesTree, 
      _modelRates, 
      families, 
      _geneTrees);
  Logger::timed << "Initial ll=" << getEvaluator().computeLikelihood() 
    << std::endl;
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
}
  
double GTSpeciesTreeOptimizer::optimizeModelRates(bool thorough)
{
  return getEvaluator().optimizeModelRates(thorough);
}

void GTSpeciesTreeOptimizer::optimize()
{
  size_t hash1 = 0;
  size_t hash2 = 0;
  unsigned int index = 0;
  _searchState.bestLL = getEvaluator().computeLikelihood();
  _searchState.farFromPlausible = true;
  /**
   *  Alternate transfer search and normal
   *  SPR search, until one does not find
   *  a better tree. Run each at least once.
   */
  rootSearch(3);
  do {
    if (index++ % 2 == 0) {
      transferSearch();
    } else {
      sprSearch(1);
    }
    if (!_searchState.farFromPlausible) {
      rootSearch(3);
    }
    hash1 = _speciesTree->getHash();
  }
  while(testAndSwap(hash1, hash2));
  rootSearch(-1);

}
double GTSpeciesTreeOptimizer::sprSearch(unsigned int radius)
{
  SpeciesSPRSearch::SPRSearch(*_speciesTree,
      getEvaluator(),
      _searchState,
      radius);
  Logger::timed << "After normal search: LL=" << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}


void GTSpeciesTreeOptimizer::onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate)
{
  getEvaluator().onSpeciesTreeChange(nodesToInvalidate);
}


void GTSpeciesTreeOptimizer::saveSpeciesTree()
{
  saveCurrentSpeciesTreeId();
}

std::string GTSpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  _speciesTree->getDatedTree().rescaleBranchLengths();
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  if (_rootLikelihoods.idToLL.size()) {
    auto newick = _speciesTree->getTree().getNewickString();
    PLLRootedTree tree(newick, false); 
    _rootLikelihoods.fillTree(tree);
    auto out = Paths::getSpeciesTreeFile(_outputDir, 
        "species_tree_llr.newick");
    tree.save(out);
  }
  return res;
}

void GTSpeciesTreeOptimizer::saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly)
{
  _speciesTree->saveToFile(str, masterRankOnly);
  if (masterRankOnly) {
    ParallelContext::barrier();
  }
}
  
double GTSpeciesTreeOptimizer::rootSearch(unsigned int maxDepth, bool thorough)
{
  _rootLikelihoods.reset();
  if (thorough) {
    _searchState.farFromPlausible = thorough;
  }
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      getEvaluator(),
      _searchState,
      maxDepth,
      &_rootLikelihoods);
  
  return _searchState.bestLL;
}
  
double GTSpeciesTreeOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
    *_speciesTree,
    getEvaluator(),
    _searchState,
    _outputDir);
  Logger::timed << "After normal search: LL=" 
    << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}
  

void GTSpeciesTreeOptimizer::reconcile(unsigned int samples)
{
  if (samples == 0) {
    return;
  }
  auto recDir = FileSystem::joinPaths(_outputDir, "reconciliations");
  FileSystem::mkdir(recDir, true);
  auto allRecDir = FileSystem::joinPaths(recDir, "all");
  FileSystem::mkdir(allRecDir, true);
  auto summariesDir = FileSystem::joinPaths(recDir, "summaries");
  FileSystem::mkdir(summariesDir, true);
  ParallelContext::barrier();
  auto &families = _geneTrees.getTrees();
  std::vector<std::string> summaryPerSpeciesEventCountsFiles;
  std::vector<std::string> summaryTransferFiles;
  for (unsigned int i = 0; i < families.size(); ++i) {
    std::vector<std::string> perSpeciesEventCountsFiles;
    std::vector<std::string> transferFiles;
    std::string geneTreesPath = FileSystem::joinPaths(allRecDir, families[i].name + std::string(".newick"));
    ParallelOfstream geneTreesOs(geneTreesPath, false);

    for (unsigned int sample = 0; sample < samples; ++ sample) {
      auto out = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_") +  std::to_string(sample) + ".xml");
      auto eventCountsFile = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_eventcount_") +  std::to_string(sample) + ".txt");
      auto perSpeciesEventCountsFile = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_perspecies_eventcount_") +  std::to_string(sample) + ".txt");
      auto transferFile = FileSystem::joinPaths(allRecDir, 
          families[i].name + std::string("_transfers_") +  std::to_string(sample) + ".txt");
      perSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
      transferFiles.push_back(transferFile);
      auto &evaluation = _evaluator->getEvaluation(i);
      evaluation.computeLogLikelihood();
      Scenario scenario;
      evaluation.inferMLScenario(scenario, true);
      scenario.saveReconciliation(out, ReconciliationFormat::RecPhyloXML, false);
      scenario.saveReconciliation(geneTreesOs, ReconciliationFormat::NewickEvents);
      scenario.saveEventsCounts(eventCountsFile, false);
      scenario.savePerSpeciesEventsCounts(perSpeciesEventCountsFile, false);
      scenario.saveTransfers(transferFile, false);
      geneTreesOs <<  "\n";
    }
    auto perSpeciesEventCountsFile = FileSystem::joinPaths(summariesDir, families[i].name + 
        std::string("_perspecies_eventcount.txt"));
    Scenario::mergePerSpeciesEventCounts(_speciesTree->getTree(),
        perSpeciesEventCountsFile,
        perSpeciesEventCountsFiles, false);
    summaryPerSpeciesEventCountsFiles.push_back(perSpeciesEventCountsFile);
    auto transferFile = FileSystem::joinPaths(summariesDir, families[i].name + 
        std::string("_transfers.txt"));
    Scenario::mergeTransfers(_speciesTree->getTree(),
        transferFile,
        transferFiles, false);
    summaryTransferFiles.push_back(transferFile);
  }
  auto totalPerSpeciesEventCountsFile = FileSystem::joinPaths(recDir, "perspecies_eventcount.txt");
  Scenario::mergePerSpeciesEventCounts(_speciesTree->getTree(),
      totalPerSpeciesEventCountsFile, 
      summaryPerSpeciesEventCountsFiles, 
      true);
  auto totalTransferFile = FileSystem::joinPaths(recDir, "transfers.txt");
  Scenario::mergeTransfers(_speciesTree->getTree(),
      totalTransferFile, 
      summaryTransferFiles, 
      true);
  ParallelContext::makeRandConsistent();
}

void GTSpeciesTreeOptimizer::optimizeDates(bool thorough)
{
  if (!_info.isDated()) {
    return; 
  }
  DatedSpeciesTreeSearch::optimizeDates(*_speciesTree,
      getEvaluator(),
      _searchState,
      getEvaluator().computeLikelihood(),
      thorough);
}



void GTSpeciesTreeOptimizer::randomizeRoot()
{
  auto &tree = _speciesTree->getDatedTree();
  unsigned int N = tree.getOrderedSpeciations().size();
  for (unsigned int i = 0; i < N; ++i) {
    auto direction = Random::getInt() % 4;
    if (SpeciesTreeOperator::canChangeRoot(*_speciesTree, direction)) {
      SpeciesTreeOperator::changeRoot(*_speciesTree, direction);
    }
  }
}


