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

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihood()
{
  return computeLikelihoodFast();
}

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihoodFast()
{
  double sumLL = 0.0;
  for (auto &evaluation: *_evaluations) {
    auto ll = evaluation->computeLogLikelihood();
    sumLL += ll;
  }
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

void GTSpeciesTreeLikelihoodEvaluator::setAlpha(double alpha)
{
  for (auto evaluation: *_evaluations) {
    evaluation->setAlpha(alpha);
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
  auto rates = *_modelRates;
  OptimizationSettings settings;
  double ll = computeLikelihood();
  Logger::timed << "[Species search] Before model rate opt, ll=" << ll << " initial rates: " << _modelRates->rates << std::endl;
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
  GTMultiEvaluatorsFunction function(*_evaluations, 
      _modelRates->info, 
      _speciesTree->getTree().getNodesNumber());
  _modelRates->rates = DTLOptimizer::optimizeParameters(
      function, 
      _modelRates->rates, 
      settings);
  ll = computeLikelihood();
  Logger::timed << "[Species search]   After model rate opt, ll=" << ll << " initial rates: " << _modelRates->rates << std::endl;
  ll = optimizeGammaRates();
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
  auto gammaCategories = _modelRates->info.gammaCategories;
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
  std::vector<double> categories(_modelRates->info.gammaCategories);
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
  auto infoCopy = _modelRates->info;
  //infoCopy.transferConstraint = TransferConstaint::NONE;
  infoCopy.originationStrategy = OriginationStrategy::UNIFORM;
  for (const auto &geneTree: _geneTrees->getTrees()) {
    auto &family = (*_families)[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    UndatedDTLMultiModel<ScaledValue> evaluation(
        speciesTree.getDatedTree(),
        mapping,
        infoCopy,
        family.startingGeneTree);
    Scenario scenario;
    // warning, this might make the random state 
    // inconsistent between the MPI ranks
    // ParallelContext::makeRandConsistent() needs to be called 
    // right after the loop
    evaluation.computeLogLikelihood();
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
  for (auto &evaluation: *_evaluations) {
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
  Logger::timed << "Initializing ccps" << std::endl;
  for (const auto &geneTree: _geneTrees.getTrees()) {
    auto &family = families[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    switch (info.model) {
    case RecModel::UndatedDL:
      _evaluations.push_back(
        std::make_shared<UndatedDLMultiModel<double> >(
        _speciesTree->getTree(),
        mapping,
        info,
        family.startingGeneTree));
      break;
    case RecModel::UndatedDTL:
      _evaluations.push_back(std::make_shared<UndatedDTLMultiModel<ScaledValue> >(
        _speciesTree->getDatedTree(),
        mapping,
        info,
        family.startingGeneTree));
      break;
    default:
      assert(false);
      break;
    }
  }
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
  Logger::timed << "Load balancing: " << averageCladesNumber / maxCladesNumber << std::endl;
  _speciesTree->addListener(this);
  ParallelContext::barrier();
  _evaluator.setEvaluations(*_speciesTree, _modelRates, families, _evaluations, _geneTrees);
  
  Logger::timed << "Initial ll=" << _evaluator.computeLikelihood() << std::endl;
 
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
}
  
double GTSpeciesTreeOptimizer::optimizeModelRates(bool thorough)
{
  return _evaluator.optimizeModelRates(thorough);
}

void GTSpeciesTreeOptimizer::optimize()
{
  size_t hash1 = 0;
  size_t hash2 = 0;
  unsigned int index = 0;
  _searchState.bestLL = _evaluator.computeLikelihood();
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
      _evaluator,
      _searchState,
      radius);
  Logger::timed << "After normal search: LL=" << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}


void GTSpeciesTreeOptimizer::onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate)
{
  for (auto &evaluation: _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
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
  _searchState.farFromPlausible = thorough;
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      _evaluator,
      _searchState,
      maxDepth,
      &_rootLikelihoods);
  
  return _searchState.bestLL;
}
  
double GTSpeciesTreeOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
    *_speciesTree,
    _evaluator,
    _searchState,
    _outputDir);
  Logger::timed << "After normal search: LL=" 
    << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}
  
void GTSpeciesTreeOptimizer::printFamilyDimensions(const std::string &outputFile)
{
  std::vector<unsigned int> localTreeSizes;
  std::vector<unsigned int> localCcpSizes;
  std::vector<unsigned int> treeSizes;
  std::vector<unsigned int> ccpSizes;
  for (auto &evaluation: _evaluations) {
    auto &ccp = evaluation->getCCP();
    localTreeSizes.push_back(ccp.getLeafNumber());
    localCcpSizes.push_back(ccp.getCladesNumber());
  }
  ParallelContext::concatenateHetherogeneousUIntVectors(
      localTreeSizes, treeSizes);
  ParallelContext::concatenateHetherogeneousUIntVectors(
      localCcpSizes, ccpSizes);
  std::vector<PairUInt> dims(treeSizes.size());
  for (unsigned int i = 0; i < treeSizes.size(); ++i) {
    dims[i].first = treeSizes[i];
    dims[i].second = ccpSizes[i];
  }
  std::sort(dims.begin(), dims.end());
  ParallelOfstream os(outputFile);
  for (auto &dim: dims) {
    os << dim.first << "," << dim.second << std::endl;
  }
}

void GTSpeciesTreeOptimizer::reconcile(unsigned int samples)
{
  if (samples == 0) {
    return;
  }
  auto recDir = FileSystem::joinPaths(_outputDir, "reconciliations");
  FileSystem::mkdir(recDir, true);
  ParallelContext::barrier();
  auto &families = _geneTrees.getTrees();

  for (unsigned int i = 0; i < families.size(); ++i) {
    std::string geneTreesPath = FileSystem::joinPaths(recDir, families[i].name + std::string(".newick"));
    ParallelOfstream geneTreesOs(geneTreesPath, false);

    for (unsigned int sample = 0; sample < samples; ++ sample) {
      auto out = FileSystem::joinPaths(recDir, families[i].name + std::string("_") +  std::to_string(sample) + ".xml");
      auto &evaluation = _evaluations[i];
      evaluation->computeLogLikelihood();
      Scenario scenario;
      evaluation->inferMLScenario(scenario, true);
      scenario.saveReconciliation(out, ReconciliationFormat::RecPhyloXML, false);
      scenario.saveReconciliation(geneTreesOs, ReconciliationFormat::NewickEvents);
      geneTreesOs <<  "\n";
    }
  }

  ParallelContext::makeRandConsistent();
}

void GTSpeciesTreeOptimizer::optimizeDates(bool thorough)
{
  if (!_info.isDated()) {
    return; 
  }
  DatedSpeciesTreeSearch::optimizeDates(*_speciesTree,
      _evaluator,
      _searchState,
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


