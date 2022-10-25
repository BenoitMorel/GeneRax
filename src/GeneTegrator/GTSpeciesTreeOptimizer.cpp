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
  Logger::timed << "[Species search] Initial rates: " << _modelRates->rates << std::endl;
  OptimizationSettings settings;
  double ll = computeLikelihood();
  Logger::timed << "[Species search] Initial likleihood " << ll << std::endl;
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
  //if (!_modelRates->info.perFamilyRates) {
    Logger::timed << "[Species search] Best rates: " << _modelRates->rates << std::endl;
  //}
    Logger::info << "New ll = " << computeLikelihood() << std::endl;
    return _modelRates->rates.getScore();
}
  
void GTSpeciesTreeLikelihoodEvaluator::getTransferInformation(PLLRootedTree &speciesTree,
    TransferFrequencies &transferFrequencies,
    PerSpeciesEvents &perSpeciesEvents)
{
  assert(false);
  /*
  // this is duplicated code from Routines...
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int labelsNumber = idToLabel.size();
  const VectorUint zeros(labelsNumber, 0);
  transferFrequencies.count = MatrixUint(labelsNumber, zeros);
  transferFrequencies.idToLabel = idToLabel;
  perSpeciesEvents = PerSpeciesEvents(speciesTree.getNodesNumber());
  auto infoCopy = _modelRates->info;
  for (const auto &geneTree: _geneTrees->getTrees()) {
    auto &family = (*_families)[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    UndatedDTLMultiModel<ScaledValue> evaluation(
        speciesTree,
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
  */
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
  Logger::info << "Load balancing: " << maxCladesNumber << " " << averageCladesNumber <<  "(" << worstFamily << ")"  << std::endl;

  Logger::timed << "Initializing ccps finished" << std::endl;
  _speciesTree->addListener(this);
  ParallelContext::barrier();
  _evaluator.setEvaluations(*_speciesTree, _modelRates, families, _evaluations, _geneTrees);
  
  Logger::timed << "Initial ll=" << _evaluator.computeLikelihood() << std::endl;
  printFamilyDimensions("temp.txt");
  
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



std::string GTSpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  return res;
}

void GTSpeciesTreeOptimizer::saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly)
{
  _speciesTree->saveToFile(str, masterRankOnly);
  if (masterRankOnly) {
    ParallelContext::barrier();
  }
}
  
double GTSpeciesTreeOptimizer::rootSearch(unsigned int maxDepth)
{
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      _evaluator,
      _searchState,
      maxDepth);
  return _searchState.bestLL;
}
  
double GTSpeciesTreeOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
    *_speciesTree,
    _evaluator,
    _searchState);
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

void GTSpeciesTreeOptimizer::perturbateDates()
{
  if (!_info.isDated()) {
    return; 
  }
  auto &tree = _speciesTree->getDatedTree();
  size_t N = tree.getOrderedSpeciations().size();
  for (size_t i = 0; i < N * 2; ++i) {
    auto rank = Random::getInt() % N;
    tree.moveUp(rank);
  }
}


void GTSpeciesTreeOptimizer::optimizeDates()
{
  if (!_info.isDated()) {
    return; 
  }
  auto &tree = _speciesTree->getDatedTree();
  auto bestLL = _evaluator.computeLikelihood();
  Logger::timed << "Optimizing dates, ll=" << bestLL << std::endl;
  optimizeDatesNaive(); 
  bestLL = _evaluator.computeLikelihood();
  Logger::timed << "Best ll so far: = " << bestLL << std::endl;
  
  for (unsigned int i = 0; i < 2; ++i) {
    bool improved = true;
    while (improved) {
      improved = false;
      auto backup = tree.getBackup();
      perturbateDates();
      auto ll = _evaluator.computeLikelihood();
      //Logger::timed << "LL after random perturbate: " << ll << std::endl;
      optimizeDatesNaive();
      ll = _evaluator.computeLikelihood();
      //Logger::timed << "LL after naive search from random perturbate: " << ll << std::endl;
      if (ll <= bestLL) {
        tree.restore(backup);
      } else {
        bestLL = ll;
        Logger::timed << "Best ll so far: = " << bestLL << std::endl;
        improved = true;
      }
    }
  } 
}


void GTSpeciesTreeOptimizer::optimizeDatesNaive()
{
  if (!_info.isDated()) {
    return; 
  }
  auto &tree = _speciesTree->getDatedTree();
  
  unsigned int max = tree.getOrderedSpeciations().size();
  auto bestLL = _evaluator.computeLikelihood();
  bool tryAgain = true;
  while (tryAgain) {
    tryAgain = false;
    Logger::timed << "new naive dated round ll= " << bestLL <<  std::endl;
    for (unsigned int rank = 0; rank < max; ++rank) {
      if (!tree.moveUp(rank)) {
        continue;
      }
      onSpeciesTreeChange(nullptr);
      auto ll = _evaluator.computeLikelihood();
      if (ll > bestLL) {
        tryAgain = true;
        bestLL = ll;
        //Logger::info << "Better ll = " << bestLL << std::endl;
        rank -= std::min((unsigned int)2, rank);
      } else {
        tree.moveUp(rank);
        onSpeciesTreeChange(nullptr);
      }
      assert(_speciesTree->getDatedTree().isConsistent());
    }
  }
    Logger::timed << "end naive dated round ll= " << bestLL <<  std::endl;
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


double GTSpeciesTreeOptimizer::plopRoot()
{
  auto ll = _evaluator.computeLikelihood();
  Logger::info << "Before root change: ll= " << ll << std::endl;
  unsigned int direction = 3;
  if (SpeciesTreeOperator::canChangeRoot(*_speciesTree, direction)) {
    SpeciesTreeOperator::changeRoot(*_speciesTree, direction);
    auto ll2 = _evaluator.computeLikelihood();
    Logger::info << "After root change: ll= " << ll2 << std::endl;
    assert(_speciesTree->getDatedTree().isConsistent());
    
  }
  if (SpeciesTreeOperator::canChangeRoot(*_speciesTree, direction)) {
    SpeciesTreeOperator::changeRoot(*_speciesTree, direction);
    auto ll2 = _evaluator.computeLikelihood();
    Logger::info << "After root change: ll= " << ll2 << std::endl;
    assert(_speciesTree->getDatedTree().isConsistent());
    
  }

  return 0.0;


}

