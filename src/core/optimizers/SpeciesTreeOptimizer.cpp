#include "SpeciesTreeOptimizer.hpp"

#include <optimizers/DTLOptimizer.hpp>
#include <IO/FileSystem.hpp>
#include <routines/Routines.hpp>
#include <algorithm>
#include <trees/TreeDuplicatesFinder.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <fstream>

static std::string getStepTag(bool fastMove)
{
  static std::string fastMoveString("[Species tree search - Fast moves]");
  static std::string slowMoveString("[Species tree search - Slow moves]");
  return fastMove ? fastMoveString : slowMoveString;
}

SpeciesTreeOptimizer::SpeciesTreeOptimizer(const std::string speciesTreeFile, 
    const Families &initialFamilies, 
    RecModel model,
    const Parameters &startingRates,
    bool perFamilyRates,
    bool userDTLRates,
    double minGeneBranchLength,
    bool pruneSpeciesTree,
    double supportThreshold,
    const std::string &outputDir,
    const std::string &execPath,
    bool fractionMissing,
    bool constrainSearch):
  _speciesTree(nullptr),
  _geneTrees(nullptr),
  _initialFamilies(initialFamilies),
  _currentFamilies(initialFamilies),
  _outputDir(outputDir),
  _execPath(execPath),
  _geneTreeIteration(1000000000), // we need to find a better way for avoiding directories collision
  _supportThreshold(supportThreshold),
  _lastRecLL(-std::numeric_limits<double>::infinity()),
  _lastLibpllLL(-std::numeric_limits<double>::infinity()),
  _bestRecLL(-std::numeric_limits<double>::infinity()),
  _bestLibpllLL(-std::numeric_limits<double>::infinity()),
  _firstOptimizeRatesCall(true),
  _userDTLRates(userDTLRates),
  _minGeneBranchLength(minGeneBranchLength),
  _pruneSpeciesTree(pruneSpeciesTree),
  _modelRates(startingRates, model, false, 1),
  _constrainSearch(constrainSearch),
  _okForClades(0),
  _koForClades(0)
{
  if (fractionMissing) {
    _fractionMissingFile = FileSystem::joinPaths(
      _outputDir,
      std::string("fractionMissing.txt"));
  }
  if (speciesTreeFile == "random") {
    _speciesTree = std::make_unique<SpeciesTree>(initialFamilies);
    setGeneTreesFromFamilies(initialFamilies);
  } else {
    _speciesTree = std::make_unique<SpeciesTree>(speciesTreeFile);
    setGeneTreesFromFamilies(initialFamilies);
  }
  _modelRates = ModelParameters(startingRates, model, perFamilyRates, _geneTrees->getTrees().size());
  //_speciesTree->saveToFile(FileSystem::joinPaths(_outputDir, "starting_species_tree.newick"), true);
  _speciesTree->addListener(this);
  std::string subsamplesPath = FileSystem::joinPaths(_outputDir, "subsamples");
  FileSystem::mkdir(FileSystem::joinPaths(_outputDir, "sub_genes_opt"), true);
  FileSystem::mkdir(subsamplesPath, true);
  saveCurrentSpeciesTreeId();
  _computeAllGeneClades();
}

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

void SpeciesTreeOptimizer::optimize(SpeciesSearchStrategy strategy,
      unsigned int sprRadius)
{
  _bestRecLL = computeRecLikelihood();
  switch (strategy) {
  case SpeciesSearchStrategy::SPR:
    for (unsigned int radius = 1; radius <= sprRadius; ++radius) {
      optimizeDTLRates();
      sprSearch(radius);
    }
    break;
  case SpeciesSearchStrategy::TRANSFERS:
    transferSearch();
    break;
  case SpeciesSearchStrategy::HYBRID:
    /**
     *  Alternate transfer search and normal
     *  SPR search, until one does not find
     *  a better tree. Run each at least once.
     */
    size_t hash1 = 0;
    size_t hash2 = 0;
    unsigned int index = 0;
    optimizeDTLRates();
    do {
      if (index++ % 2 == 0) {
        transferSearch();
      } else {
        sprSearch(1);
      }
      hash1 = _speciesTree->getHash();
    }
    while(testAndSwap(hash1, hash2));
    rootExhaustiveSearch();
    break;
  }
}


SpeciesTreeOptimizer::~SpeciesTreeOptimizer()
{
  _speciesTree->removeListener(this);
}
  
void SpeciesTreeOptimizer::rootExhaustiveSearchAux(SpeciesTree &speciesTree, 
    PerCoreGeneTrees &geneTrees, 
    RecModel model, 
    bool doOptimizeGeneTrees, 
    std::vector<unsigned int> &movesHistory, 
    std::vector<unsigned int> &bestMovesHistory, 
    double &bestLL, 
    unsigned int &visits)
{
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  for (auto direction: moves) {
    if (SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      movesHistory.push_back(direction);
      SpeciesTreeOperator::changeRoot(speciesTree, direction);
      unsigned int geneRadius = doOptimizeGeneTrees ? 1 : 0;
      double ll = computeLikelihood(geneRadius);
      visits++;
      if (ll > bestLL) {
        bestLL = ll;
        bestMovesHistory = movesHistory; 
        Logger::info << "Found better root " << ll << std::endl;
      }
      rootExhaustiveSearchAux(speciesTree, 
          geneTrees, 
          model, 
          doOptimizeGeneTrees, 
          movesHistory, 
          bestMovesHistory, 
          bestLL, 
          visits);
      SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
      movesHistory.pop_back();
    }
  }
}

double SpeciesTreeOptimizer::rootExhaustiveSearch(bool doOptimizeGeneTrees)
{
  Logger::timed << "Root exhaustive search " << std::endl;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  unsigned int geneRadius = doOptimizeGeneTrees ? 1 : 0;
  double bestLL = computeLikelihood(geneRadius);
  unsigned int visits = 1;
  movesHistory.push_back(0);
  rootExhaustiveSearchAux(*_speciesTree, 
      *_geneTrees, 
      _modelRates.model, 
      doOptimizeGeneTrees, 
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits); 
  movesHistory[0] = 1;
  rootExhaustiveSearchAux(*_speciesTree, 
      *_geneTrees, 
      _modelRates.model, 
      doOptimizeGeneTrees, 
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits); 
  assert (visits == 2 * _speciesTree->getTree().getLeavesNumber() - 3);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(*_speciesTree, bestMovesHistory[i]);
  }
  return bestLL;
}



bool SpeciesTreeOptimizer::testPruning(unsigned int prune,
    unsigned int regraft,
    double refApproxLL, 
    unsigned int hash1)
{
  auto wrongClades1 = _unsupportedCladesNumber();
  bool tryApproxFirst = Enums::implementsApproxLikelihood(_modelRates.model);
  bool check = false;
  // Apply the move
  auto rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
  _stats.testedTrees++;
  bool canTestMove = true;
  // Discard bad moves with an approximation of the likelihood function
  double approxRecLL;
  bool needFullRollback = false;
  auto wrongClades2 = _unsupportedCladesNumber();
  
  if (_constrainSearch && wrongClades2 > wrongClades1) {
    canTestMove = false;
    _koForClades++;
  } else {
    _okForClades++;
  }
  if (tryApproxFirst && canTestMove) {
    approxRecLL = computeApproxRecLikelihood();
    if (approxRecLL - _bestRecLL < 0.0) {
      canTestMove = false;
    } else {
      needFullRollback = true;
    }
  }
  if (canTestMove) {
    // we really test the move
    _lastRecLL = computeRecLikelihood();
    if (_lastRecLL > _bestRecLL) {
      // Better tree found! keep it and return
      newBestTreeCallback();
      return true;
    }
  }
  // we do not keep the tree
  SpeciesTreeOperator::reverseSPRMove(*_speciesTree, prune, rollback);
  if (needFullRollback) {
    for (auto &evaluation: _evaluations) {
      evaluation->rollbackToLastState();
    }
  }
  // ensure that we correctly reverted
  if (check) {
    auto hash2 = _speciesTree->getNodeIndexHash(); 
    assert(hash1 == hash2);
    if (tryApproxFirst && !canTestMove) {
      auto approxRevertedLL = computeApproxRecLikelihood();
      assert(fabs(refApproxLL - approxRevertedLL) < 0.1);
    } else {
      auto revertedLL = computeRecLikelihood();
      assert(fabs(revertedLL - _bestRecLL) < 0.1);
    }
  }
  return false;
}

struct TransferMove {
  unsigned int prune;
  unsigned int regraft;
  double transfers;
  TransferMove(): prune(0), regraft(0), transfers(0.0) {
  }
  TransferMove(unsigned int p, unsigned int r, double t): prune(p), regraft(r), transfers(t) {
  }
  bool operator < (const TransferMove& tm) const
  {
    if (transfers != tm.transfers) {
      return transfers > tm.transfers;
    } else if (regraft != tm.regraft) {
      return regraft > tm.regraft;
    } else {
      return prune > tm.prune;
    }
  }
  bool operator ==(const TransferMove& obj) const
  {
    return (obj.prune == prune) && (obj.regraft == regraft) && (obj.transfers == transfers); 
  }
};

static unsigned int hashints(unsigned int a, unsigned int b)
{
  return (a + b) * (a + b + 1) / 2  + b;
}

namespace std {
template<>
struct hash<TransferMove>
{
  size_t
    operator()(const TransferMove & obj) const
    {
      return hash<int>()(static_cast<int>(
            hashints(hashints(obj.prune, obj.regraft), obj.transfers)));
    }
};
}



struct MovesBlackList {
  std::unordered_set<TransferMove> _blacklist;
  bool isBlackListed(const TransferMove &move) { return _blacklist.find(move) != _blacklist.end();}
  void blacklist(const TransferMove &move) { _blacklist.insert(move); }
};

static std::unordered_map<std::string, double> getCoverage(const std::string &path)
{
  std::ifstream is(path);
  std::unordered_map<std::string, double> res;
  std::string line;
  std::getline(is, line);
  while (std::getline(is, line))
  {
    std::istringstream iss(line);
    std::string species;
    double coverage;
    iss >> species >> coverage;
    species.pop_back();
    res.insert({species, coverage});
  }
  return res;
}

double SpeciesTreeOptimizer::fastTransfersRound(MovesBlackList &blacklist)
{
  unsigned int minTransfers = 1;
  _bestRecLL = computeRecLikelihood();
  auto hash1 = _speciesTree->getNodeIndexHash(); 
  auto refApproxLL = computeApproxRecLikelihood();
  TransferFrequencies frequencies;
  std::string speciesTreeFile(FileSystem::joinPaths(_outputDir, "speciesTreeTemp.newick"));
  saveCurrentSpeciesTreePath(speciesTreeFile, true);
  ParallelContext::barrier();
  Logger::timed << "Start inferring transfers..." << std::endl;
  Routines::getTransfersFrequencies(speciesTreeFile,
    _currentFamilies,
    _modelRates,
    _pruneSpeciesTree,
    frequencies,
    _outputDir);
  auto coverages = getCoverage(_outputDir + "/perSpeciesCoverage.txt");
  unsigned int transfers = 0;
  ParallelContext::barrier();
  std::unordered_map<std::string, unsigned int> labelsToIds;
  _speciesTree->getLabelsToId(labelsToIds);
  std::vector<TransferMove> transferMoves;
  static int plop = 0;
  plop++;
  for (auto entry: frequencies) {
    transfers += entry.second;
    if (entry.second >= minTransfers) {
      std::string key1, key2;
      Routines::getLabelsFromTransferKey(entry.first, key1, key2);
      unsigned int prune = labelsToIds[key2];
      unsigned int regraft = labelsToIds[key1];
      // HERE
      if (SpeciesTreeOperator::canApplySPRMove(*_speciesTree, prune, regraft)) {
        TransferMove move(prune, regraft, entry.second);
        double factor = 1.0;
        if (plop % 2 == 0) {
          if (coverages.find(key1) != coverages.end()) {
            assert(coverages[key1] != 0.0);
            factor /= coverages[key1];
          }
          if (coverages.find(key2) != coverages.end()) {
            assert(coverages[key2] != 0.0);
            factor /= coverages[key2];
          }
        }
        if (!blacklist.isBlackListed(move)) {
          transferMoves.push_back(TransferMove(prune, regraft, factor * entry.second)); 
        }
      }
    }
  }
  Logger::timed << "  Inferred transfers:" << transfers << std::endl;
  std::sort(transferMoves.begin(), transferMoves.end());
  unsigned int index = 0;
  const unsigned int stopAfterFailures = 50;
  const unsigned int stopAfterImprovements = 50;
  const unsigned int minTrial = _speciesTree->getTree().getNodesNumber();
  unsigned int failures = 0;
  unsigned int improvements = 0;
  for (auto &transferMove: transferMoves) {
    index++;
    _stats.testedTransfers++;
    if (SpeciesTreeOperator::canApplySPRMove(*_speciesTree, transferMove.prune, transferMove.regraft)) {
      blacklist.blacklist(transferMove);
      if (testPruning(transferMove.prune, transferMove.regraft, refApproxLL, hash1)) {
        _stats.acceptedTransfers++;
        failures = 0;
        improvements++;
        Logger::info << "  better tree (transfers:" << transferMove.transfers << ", trial: " << index << ", ll=" << _bestRecLL << ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() << ")"   << std::endl;
        // we enough improvements to recompute the new transfers
        hash1 = _speciesTree->getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        refApproxLL = computeApproxRecLikelihood();
      } else {
        failures++;
      }
      bool stop = index > minTrial && failures > stopAfterFailures;
      stop |= improvements > stopAfterImprovements;
      if (stop) {
        return _bestRecLL;
      }
    }  
  }
  return _bestRecLL;
}

double SpeciesTreeOptimizer::fastSPRRound(unsigned int radius)
{
  _bestRecLL = computeRecLikelihood();
  auto hash1 = _speciesTree->getNodeIndexHash(); 
  auto refApproxLL = computeApproxRecLikelihood();

  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  //assert (fabs(_bestRecLL - refApproxLL) < 0.01);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      if (testPruning(prune, regraft, refApproxLL, hash1)) {
        Logger::timed << "\tbetter tree (LL=" 
          << _bestRecLL << ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() << ")"<< std::endl;
        hash1 = _speciesTree->getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        refApproxLL = computeApproxRecLikelihood();
      }
    }
  }
  return _bestRecLL;
}


struct less_than_evaluatedmove
{
  inline bool operator() (const EvaluatedMove& e1, const EvaluatedMove& e2)
  {
    return e1.ll > e2.ll;
  }
};

std::vector<EvaluatedMove> SpeciesTreeOptimizer::getSortedCandidateMoves(unsigned int speciesRadius) 
{
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  std::vector<EvaluatedMove> evaluatedMoves;
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, speciesRadius, regrafts);
    for (auto regraft: regrafts) {
      unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
      EvaluatedMove em;
      em.prune = prune;
      em.regraft = regraft;
      em.ll = computeRecLikelihood();
      evaluatedMoves.push_back(em);
      SpeciesTreeOperator::reverseSPRMove(*_speciesTree, em.prune, rollback);
    }
  }
  std::sort(evaluatedMoves.begin(), evaluatedMoves.end(), less_than_evaluatedmove());
  return evaluatedMoves;
}
  
struct ReferenceLikelihood {
  unsigned int radius;
  double refLikelihood;
  double tolerance;
};

double SpeciesTreeOptimizer::slowSPRRound(unsigned int speciesRadius, double bestLL)
{
  const double jointLikelihoodEpsilon = -10;
  const unsigned int maxMovesToTry = 20;
  const unsigned int maxGeneRadius = 1;

  Logger::timed << getStepTag(false) 
    << " Starting new SPR round from tree hash=" << _speciesTree->getHash() << std::endl;
  std::vector<ReferenceLikelihood> referenceLikelihoods;
  for (unsigned int currentRadius = 1; currentRadius <= maxGeneRadius; ++currentRadius) {
    ReferenceLikelihood ref;
    ref.radius = currentRadius;
    ref.refLikelihood = computeLikelihood(currentRadius);
    ref.tolerance = (currentRadius == maxGeneRadius ? 0 : jointLikelihoodEpsilon);
    referenceLikelihoods.push_back(ref); 
  }
  Logger::timed << getStepTag(false) << "   Slow round from tree hash=" << _speciesTree->getHash()
    << " joint ll= " <<  referenceLikelihoods.back().refLikelihood << std::endl;
  auto sortedCandidateMoves = getSortedCandidateMoves(speciesRadius);
  unsigned int movesToTry = std::min(maxMovesToTry, static_cast<unsigned int>(sortedCandidateMoves.size()));
  for (unsigned int i = 0; i < movesToTry; ++i) {
    auto &em = sortedCandidateMoves[i];
    unsigned int rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, em.prune, em.regraft);
    bool isBetter = true;
    double newBestLL = -std::numeric_limits<double>::infinity();
    assert(referenceLikelihoods.size());
    for (auto &ref: referenceLikelihoods) {
      newBestLL = computeLikelihood(ref.radius);
      if (newBestLL < ref.refLikelihood + ref.tolerance) {
        isBetter = false;
        break;
      } 
    }
    if (isBetter) {
      Logger::timed << getStepTag(false) << "   Found better tree hash=" << _speciesTree->getHash() 
        << " ll=" << newBestLL << " (previous ll = " << referenceLikelihoods.back().refLikelihood << ")" << std::endl;
      newBestTreeCallback();
      return newBestLL;
    }
    SpeciesTreeOperator::reverseSPRMove(*_speciesTree, em.prune, rollback);
  } 
  return bestLL;
}

double SpeciesTreeOptimizer::transferSearch()
{
  _stats.reset();
  auto bestLL = computeRecLikelihood();
  Logger::timed << getStepTag(true) << " Starting species transfer search, bestLL=" 
    << bestLL <<  ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() 
    << ")" <<std::endl;
  double newLL = bestLL;
  MovesBlackList blacklist;
  int index = 0;
  do {
    if (index > 0) {
      bestLL = optimizeDTLRates();
    }
    newLL = fastTransfersRound(blacklist);
    Logger::info << "  Accepted: " << _okForClades << std::endl;
    Logger::info << "  Rejected: " << _koForClades << std::endl;
    index++;
  } while (newLL - bestLL > 1.0);
  if (index == 1) {
      bestLL = optimizeDTLRates();
  }
  Logger::timed << "After transfer search: " << newLL << std::endl;
  Logger::info << _stats << std::endl; 
  saveCurrentSpeciesTreeId();
  _stats.reset();
  return newLL;
}

double SpeciesTreeOptimizer::sprSearch(unsigned int radius, bool doOptimizeGeneTrees)
{
  _stats.reset();
  unsigned int geneRadius = doOptimizeGeneTrees ? 1 : 0;
  double bestLL = doOptimizeGeneTrees ? computeLikelihood(geneRadius) : computeRecLikelihood();
  Logger::timed << getStepTag(!doOptimizeGeneTrees) << " Starting species SPR search, radius=" 
    << radius << ", bestLL=" << bestLL  << ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() << ")" <<std::endl;
  double newLL = bestLL;
  do {
    bestLL = newLL;
    if (doOptimizeGeneTrees) {
      newLL = slowSPRRound(radius, bestLL); 
    } else {
      newLL = fastSPRRound(radius);
    }
  } while (newLL - bestLL > 0.001);
  Logger::timed << "After normal search: " << bestLL << std::endl;
  Logger::info << _stats << std::endl;
  saveCurrentSpeciesTreeId();
  Logger::info << "  Accepted: " << _okForClades << std::endl;
  Logger::info << "  Rejected: " << _koForClades << std::endl;
  return newLL;
}
  
ModelParameters SpeciesTreeOptimizer::computeOptimizedRates() 
{
  if (_userDTLRates) {
    return _modelRates;
  }
  auto rates = _modelRates;
  OptimizationSettings settings;
  settings.lineSearchMinImprovement = 10.0;
  settings.minAlpha = 0.01;
  double ll = computeRecLikelihood();
  settings.optimizationMinImprovement = std::max(3.0, ll / 1000.0);
  rates =  DTLOptimizer::optimizeModelParameters(_evaluations, !_firstOptimizeRatesCall, rates, settings);
  _firstOptimizeRatesCall = true;
  return rates;
}
  
double SpeciesTreeOptimizer::optimizeDTLRates()
{
  if (_userDTLRates) {
    return computeRecLikelihood();
  }
  Logger::timed << "optimize rates " << std::endl;
  _modelRates = computeOptimizedRates();
  unsigned int i = 0;
  for (auto &evaluation: _evaluations) {
    evaluation->setRates(_modelRates.getRates(i++));
  }
  Logger::timed << "optimize rates done (LL=" << computeRecLikelihood() << ")" << std::endl;
  if (!_modelRates.perFamilyRates) {
    Logger::timed << " Best rates: " << _modelRates.rates << std::endl;
  }
  return computeRecLikelihood();
}
  
std::string SpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = FileSystem::joinPaths(_outputDir, name);
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  return res;
}

void SpeciesTreeOptimizer::saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly)
{
  _speciesTree->saveToFile(str, masterRankOnly);
}

double SpeciesTreeOptimizer::optimizeGeneTrees(unsigned int radius)
{
  saveCurrentSpeciesTreeId("proposal_species_tree.newick");
  std::string speciesTree = FileSystem::joinPaths(_outputDir, "proposal_species_tree.newick");
  auto recOpt = RecOpt::Simplex;
  bool rootedGeneTree = true;
  double recWeight = 1.0;
  bool useSplitImplem = true;
  long int sumElapsedSPR = 0;
  auto rates = _modelRates;
  std::string resultName = "proposals";
  unsigned int iterationsNumber = 1;
  bool inPlace = false; 
  bool perFamilyDTLRates = false;
  bool madRooting = false;
  assert(perFamilyDTLRates == false);
  if (radius == 1) {
    iterationsNumber = 2;
  }
  for (unsigned i = 0; i < iterationsNumber; ++i) {
    Logger::mute();
    Routines::optimizeGeneTrees(_currentFamilies, 
      _modelRates.model, rates.rates, _outputDir, resultName, 
      _execPath, speciesTree, recOpt, perFamilyDTLRates, rootedGeneTree, 
      madRooting,
      _supportThreshold, recWeight, true, true, radius, _geneTreeIteration, 
        useSplitImplem, sumElapsedSPR, inPlace);
    _geneTreeIteration++;
    Logger::unmute();
    setGeneTreesFromFamilies(_currentFamilies);
    if (i < iterationsNumber - 1) {
      rates = computeOptimizedRates();
    }
  }
  Routines::gatherLikelihoods(_currentFamilies, _lastLibpllLL, _lastRecLL);
  return _lastLibpllLL + _lastRecLL;
}
 
void SpeciesTreeOptimizer::revertGeneTreeOptimization()
{
  _currentFamilies = _initialFamilies;
  setGeneTreesFromFamilies(_currentFamilies);
}

double SpeciesTreeOptimizer::computeLikelihood(unsigned int geneSPRRadius)
{
  if (geneSPRRadius >= 1) {
    double res = optimizeGeneTrees(geneSPRRadius);
    revertGeneTreeOptimization();
    return res;
  } else {
    _lastRecLL = computeRecLikelihood();
    return _lastRecLL;
  }
}

double SpeciesTreeOptimizer::computeRecLikelihood()
{
  double ll = 0.0;
  for (auto &evaluation: _evaluations) {
    ll += evaluation->evaluate(false);
  }
  ParallelContext::sumDouble(ll);
  _stats.exactLikelihoodCalls++;
  return ll;
}

double SpeciesTreeOptimizer::computeApproxRecLikelihood()
{
  double ll = 0.0;
  for (auto &evaluation: _evaluations) {
    ll += evaluation->evaluate(true);
  }
  ParallelContext::sumDouble(ll);
  _stats.approxLikelihoodCalls++;
  return ll;
}

void SpeciesTreeOptimizer::newBestTreeCallback()
{
  saveCurrentSpeciesTreeId();
  _stats.acceptedTrees++;
  _bestLibpllLL = _lastLibpllLL;
  _bestRecLL = _lastRecLL;
}
  
void SpeciesTreeOptimizer::setGeneTreesFromFamilies(const Families &families)
{
  _geneTrees = std::make_unique<PerCoreGeneTrees>(families, true);
  updateEvaluations();
}
  
void SpeciesTreeOptimizer::updateEvaluations()
{
  assert(_geneTrees);
  auto &trees = _geneTrees->getTrees();
  _evaluations.resize(trees.size());
  for (unsigned int i = 0; i < trees.size(); ++i) {
    auto &tree = trees[i];
    _evaluations[i] = std::make_shared<ReconciliationEvaluation>(_speciesTree->getTree(), *tree.geneTree, tree.mapping, _modelRates.model, false, _minGeneBranchLength, _pruneSpeciesTree, _fractionMissingFile);
    _evaluations[i]->setRates(_modelRates.getRates(i));
    _evaluations[i]->setPartialLikelihoodMode(PartialLikelihoodMode::PartialSpecies);
    //_evaluations[i]->enableMADRooting(true);
  }
}

void SpeciesTreeOptimizer::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  for (auto &evaluation: _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}

void SpeciesTreeOptimizer::likelihoodsSnapshot()
{
  std::string filename = _outputDir + "/snapshot_core" + std::to_string(ParallelContext::getRank()); 
  std::ofstream os(filename);
  auto &trees = _geneTrees->getTrees();
  for (unsigned int i = 0; i < trees.size(); ++i) {
    os << trees[i].name << " ";
    os << _evaluations[i]->evaluate() << std::endl;
  }
}

std::string getCladesSetPath(const std::string &outputDir,
    int rank)
{
  std::string basePath = "clades_" + std::to_string(rank) + ".txt";
  return FileSystem::joinPaths(outputDir, basePath);
}

void SpeciesTreeOptimizer::_computeAllGeneClades()
{
  Logger::timed << "Computing gene clades..." << std::endl;
  ParallelContext::barrier();
  
  // Compute local clades
  auto speciesLabelToInt = _speciesTree->getTree().getLabelToIntMap();
  CladeSet allClades; 
  for (auto &tree: _geneTrees->getTrees()) {
    auto cladesSet = Clade::buildCladeSet(*tree.geneTree,
        tree.mapping,
        speciesLabelToInt);
    allClades.insert(cladesSet.begin(), cladesSet.end());
  }
  // Write local clades
  std::ofstream os(
      getCladesSetPath(_outputDir, ParallelContext::getRank()));
  for (auto clade: allClades) {
    os << clade << " ";
  }
  os.close();
  ParallelContext::barrier();
  // Load all clades
  _geneClades.clear();
  for (unsigned int rank = 0; rank < ParallelContext::getSize(); ++rank) {
    std::ifstream is(getCladesSetPath(_outputDir, rank));
    unsigned int clade = 0;
    while (is >> clade) {
      _geneClades.insert(clade);
    }
  }
  assert(ParallelContext::isIntEqual(_geneClades.size()));
  Logger::timed << "Number number of supported bipartitions: " << _geneClades.size()  << std::endl;
}


unsigned int SpeciesTreeOptimizer::_unsupportedCladesNumber()
{
  auto speciesClades = Clade::buildCladeSet(_speciesTree->getTree());
  unsigned int intersectionSize =  intersection_size(speciesClades, 
      _geneClades);
  return speciesClades.size() - intersectionSize;
}

