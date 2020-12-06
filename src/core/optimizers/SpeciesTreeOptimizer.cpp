#include "SpeciesTreeOptimizer.hpp"

#include <util/Paths.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <IO/FileSystem.hpp>
#include <routines/Routines.hpp>
#include <algorithm>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <fstream>
#include <NJ/MiniNJ.hpp>
#include <NJ/NeighborJoining.hpp>
#include <cstdio>
#include <support/ICCalculator.hpp>

SpeciesTreeOptimizer::SpeciesTreeOptimizer(const std::string speciesTreeFile, 
    const Families &initialFamilies, 
    const RecModelInfo &recModelInfo,
    const Parameters &startingRates,
    bool userDTLRates,
    const std::string &outputDir,
    bool constrainSearch):
  _speciesTree(nullptr),
  _geneTrees(nullptr),
  _initialFamilies(initialFamilies),
  _outputDir(outputDir),
  _lastRecLL(-std::numeric_limits<double>::infinity()),
  _bestRecLL(-std::numeric_limits<double>::infinity()),
  _firstOptimizeRatesCall(true),
  _userDTLRates(userDTLRates),
  _modelRates(startingRates, 1, recModelInfo),
  _constrainSearch(constrainSearch),
  _okForClades(0),
  _koForClades(0),
  _hardToFindBetter(false),
  _optimizationCriteria(ReconciliationLikelihood)
{

  _modelRates.info.perFamilyRates = false; // we set it back a few
                                           // lines later
  if (speciesTreeFile == "random") {
    _speciesTree = std::make_unique<SpeciesTree>(initialFamilies);
    setGeneTreesFromFamilies(initialFamilies);
  } else {
    _speciesTree = std::make_unique<SpeciesTree>(speciesTreeFile);
    setGeneTreesFromFamilies(initialFamilies);
  }
  _modelRates = ModelParameters(startingRates, 
      _geneTrees->getTrees().size(),
      recModelInfo);
  _speciesTree->addListener(this);
  saveCurrentSpeciesTreeId();
  _computeAllGeneClades();
  _hardToFindBetter |= _unsupportedCladesNumber() <= std::max<unsigned int>(
      _speciesTree->getTree().getNodesNumber() / 10, 1);
  if (_hardToFindBetter) {
    Logger::info << "Hard-to-find-better mode" << std::endl;
  }
}

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

void SpeciesTreeOptimizer::optimize(SpeciesSearchStrategy strategy,
      unsigned int sprRadius,
      OptimizationCriteria criteria)
{
  setOptimizationCriteria(criteria);
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
    rootSearch(3);
    transferSearch();
    rootSearch();
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
    if (_hardToFindBetter) {
      optimizeDTLRates();
      _bestRecLL = computeRecLikelihood();
      rootSearch(3);
    }
    do {
      if (index++ % 2 == 0) {
        transferSearch();
      } else {
        sprSearch(sprRadius);
      }
      if (_hardToFindBetter) {
        rootSearch(3);
      }
      hash1 = _speciesTree->getHash();
    }
    while(testAndSwap(hash1, hash2));
    rootSearch(5);
    break;
  }
  setOptimizationCriteria(ReconciliationLikelihood);
}


SpeciesTreeOptimizer::~SpeciesTreeOptimizer()
{
  _speciesTree->removeListener(this);
}
  
void SpeciesTreeOptimizer::rootSearchAux(SpeciesTree &speciesTree, 
    PerCoreGeneTrees &geneTrees, 
    RecModel model, 
    std::vector<unsigned int> &movesHistory, 
    std::vector<unsigned int> &bestMovesHistory, 
    double &bestLL, 
    unsigned int &visits,
    unsigned int maxDepth)
{
  if (movesHistory.size() > maxDepth) {
    return;
  }
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  for (auto direction: moves) {
    if (SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      movesHistory.push_back(direction);
      beforeTestCallback(); 
      SpeciesTreeOperator::changeRoot(speciesTree, direction);
      optimizeGeneRoots();
      double ll = computeRecLikelihood();
      visits++;
      unsigned int plop = 0;
      if (ll > bestLL) {
        bestLL = ll;
        bestMovesHistory = movesHistory; 
        Logger::info << "Found better root " << ll << std::endl;
        plop = 3;
      }
      rootSearchAux(speciesTree, 
          geneTrees, 
          model, 
          movesHistory, 
          bestMovesHistory, 
          bestLL, 
          visits,
          maxDepth + plop);
      SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
      rollbackCallback();
      movesHistory.pop_back();
    }
  }
}

double SpeciesTreeOptimizer::rootSearch(unsigned int maxDepth)
{
  Logger::timed << "Root search with depth=" << maxDepth << std::endl;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  double bestLL = computeRecLikelihood();
  unsigned int visits = 1;
  movesHistory.push_back(0);
  rootSearchAux(*_speciesTree, 
      *_geneTrees, 
      _modelRates.info.model, 
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits, 
      maxDepth); 
  movesHistory[0] = 1;
  rootSearchAux(*_speciesTree, 
      *_geneTrees, 
      _modelRates.info.model, 
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits,
      maxDepth); 
  //assert (visits == 2 * _speciesTree->getTree().getLeavesNumber() - 3);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(*_speciesTree, bestMovesHistory[i]);
  }
  optimizeGeneRoots();
  return bestLL;
}



bool SpeciesTreeOptimizer::testPruning(unsigned int prune,
    unsigned int regraft)
{
  beforeTestCallback();
  // Apply the move
  auto rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
  bool canTestMove = true;
  const bool geneRootOpt = _modelRates.info.rootedGeneTree; 
  double geneRootApproxLL = 0.0;
  if (canTestMove && geneRootOpt && _averageGeneRootDiff.isSignificant()) {
    // approximate the likelihood with none exhaustive
    // ML gene tree root search. Decide whether we can already
    // discard the move
    geneRootApproxLL = computeRecLikelihood();
    auto epsilon = 2.0 * _averageGeneRootDiff.getAverage();
    canTestMove &= (geneRootApproxLL + epsilon > _bestRecLL );
  }
  if (!_averageGeneRootDiff.isSignificant()) {
    geneRootApproxLL = computeRecLikelihood();
  }
  if (canTestMove) {
    // we test the move with exact likelihood
    _okForClades++;
    optimizeGeneRoots();
    _lastRecLL = computeRecLikelihood();
    _averageGeneRootDiff.addValue(_lastRecLL - geneRootApproxLL);
    if (_lastRecLL > _bestRecLL) {
      // Better tree found! keep it and return
      // without rollbacking
      newBestTreeCallback();
      return true;
    }
  } else {
    _koForClades++;
  }
  // the tree is not better, rollback
  SpeciesTreeOperator::reverseSPRMove(*_speciesTree, prune, rollback);
  rollbackCallback();
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


double SpeciesTreeOptimizer::transferRound(MovesBlackList &blacklist,
    bool &maxImprovementsReached)
{
  maxImprovementsReached = false;
  unsigned int reconciliationSamples = 0;
  unsigned int minTransfers = 1;
  _bestRecLL = computeRecLikelihood();
  auto hash1 = _speciesTree->getNodeIndexHash(); 
  TransferFrequencies frequencies;
  auto  speciesTreeFile = 
    Paths::getSpeciesTreeFile(_outputDir, "speciesTreeTemp.newick");
  saveCurrentSpeciesTreePath(speciesTreeFile, true);
  ParallelContext::barrier();
  Logger::timed << "Start inferring transfers..." << std::endl;
  Routines::getTransfersFrequencies(speciesTreeFile,
    _initialFamilies,
    _modelRates,
    reconciliationSamples,
    frequencies);
  PerSpeciesEvents perSpeciesEvents;
  const bool forceTransfers = true;
  Routines::getPerSpeciesEvents(speciesTreeFile,
    _initialFamilies,
    _modelRates,
    reconciliationSamples,
    perSpeciesEvents,
    forceTransfers);
  unsigned int speciesNumber = _speciesTree->getTree().getNodesNumber();
  std::vector<double> speciesFrequencies;
  for (unsigned int e = 0; e < speciesNumber; ++e) {
    auto &speciesEvents = perSpeciesEvents.events[e];
    speciesFrequencies.push_back(speciesEvents.speciesFrequency());
  }
  Logger::timed << "Start generating transfer SPR moves......" << std::endl;
  unsigned int transfers = 0;
  ParallelContext::barrier();
  std::unordered_map<std::string, unsigned int> labelsToIds;
  _speciesTree->getLabelsToId(labelsToIds);
  std::vector<TransferMove> transferMoves;
  for (unsigned int from = 0; from < frequencies.count.size(); ++from) {
    for (unsigned int to = 0; to < frequencies.count.size(); ++to) {
      auto regraft = labelsToIds[frequencies.idToLabel[from]];
      auto prune = labelsToIds[frequencies.idToLabel[to]];
      auto count = frequencies.count[from][to];
      transfers += count;
      if (count < minTransfers) {
        continue;
      }
      if (SpeciesTreeOperator::canApplySPRMove(*_speciesTree, prune, regraft)) {
        TransferMove move(prune, regraft, count);
        double factor = 1.0;
        if (_modelRates.info.pruneSpeciesTree) {
          factor /= (1.0 + sqrt(speciesFrequencies[prune]));
          factor /= (1.0 + sqrt(speciesFrequencies[regraft]));
        }
        if (!blacklist.isBlackListed(move)) {
          transferMoves.push_back(TransferMove(prune, regraft, factor * count)); 
        }
      }
    }
  }
  Logger::timed << "  Inferred transfers:" << transfers << std::endl;
  std::sort(transferMoves.begin(), transferMoves.end());
  unsigned int index = 0;
  const unsigned int stopAfterFailures = 50u;
  const unsigned int stopAfterImprovements = std::max(15u, speciesNumber / 4);
  const unsigned int minTrial = std::max(50u, speciesNumber / 2);
  unsigned int failures = 0;
  unsigned int improvements = 0;
  std::unordered_set<unsigned int> alreadyPruned;
  Logger::timed << "Start the search..." << std::endl;
  unsigned int trials = 0;
  for (auto &transferMove: transferMoves) {
    index++;
    if (alreadyPruned.find(transferMove.prune) != alreadyPruned.end()) {
      continue;
    }
    if (SpeciesTreeOperator::canApplySPRMove(*_speciesTree, transferMove.prune, transferMove.regraft)) {
      blacklist.blacklist(transferMove);
      trials++;
      if (testPruning(transferMove.prune, transferMove.regraft)) {
        failures = 0;
        improvements++;
        alreadyPruned.insert(transferMove.prune);
        Logger::info << "  better tree (transfers:" << transferMove.transfers << ", trial: " << trials << ", ll=" << _bestRecLL << ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() << ")"   << std::endl;
        // we enough improvements to recompute the new transfers
        hash1 = _speciesTree->getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        if (_hardToFindBetter) {
          veryLocalSearch(transferMove.prune);
        }
      } else {
        failures++;
      }
      bool stop = index > minTrial && failures > stopAfterFailures;
      maxImprovementsReached = improvements > stopAfterImprovements;
      stop |= maxImprovementsReached;
      if (stop) {
        if (!_hardToFindBetter && !maxImprovementsReached) {
          Logger::timed << "Switch to hardToFindBetter mode" << std::endl;
          _hardToFindBetter = true;
        }
        return _bestRecLL;
      }
    }  
  }
  return _bestRecLL;
}

std::vector<double> SpeciesTreeOptimizer::_getSupport()
{
  std::string temp = FileSystem::joinPaths(_outputDir, "tmp");
  std::vector<double> idToSupport;
  ICCalculator::computeScores(_speciesTree->getTree(),
      _initialFamilies,
      true,
      temp,
      idToSupport);
  for (auto node: _speciesTree->getTree().getLeaves()) {
    idToSupport[node->node_index] = 1.0;
  }
  return idToSupport;
}

double SpeciesTreeOptimizer::fastSPRRound(unsigned int radius)
{
  Logger::timed << "Before SPR Round radius=" << radius << std::endl;
  Logger::timed << "Start SPR Round radius=" << radius << std::endl;
  _bestRecLL = computeRecLikelihood();
  auto hash1 = _speciesTree->getNodeIndexHash(); 

  auto supportValues = std::vector<double>();//_getSupport();
  double maxSupport = 0.2; // ignored for now
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, 
      prunes,
      supportValues,
      maxSupport);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      /*
      bool test = false;
      test |= supportValues[prune] < maxSupport;
      test |= supportValues[_speciesTree->getNode(prune)->parent->node_index] < maxSupport;
      test |= supportValues[regraft] < maxSupport;
      */
      if (testPruning(prune, regraft)) {
        Logger::timed << "\tbetter tree (LL=" 
          << _bestRecLL << ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() << ")"<< std::endl;
        auto pruneNode = _speciesTree->getNode(prune);
        
        Logger::info << pruneNode->label 
          << " " << _speciesTree->getNode(regraft)->label 
        /*
          << " " << supportValues[prune]  
          << " " << supportValues[pruneNode->parent->node_index]; 
          if (pruneNode->left) {
            Logger::info <<  " (" 
              << supportValues[pruneNode->left->node_index] 
              << ","
              << supportValues[pruneNode->right->node_index] 
              << ")";
          }
          Logger::info << " " << supportValues[regraft] 
          */
          << std::endl;
        hash1 = _speciesTree->getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        _bestRecLL = veryLocalSearch(prune);
      }
    }
  }
  return _bestRecLL;
}

double SpeciesTreeOptimizer::veryLocalSearch(unsigned int spid)
{
  const unsigned int radius = 2;
  std::vector<unsigned int> prunes;
  prunes.push_back(spid);
  unsigned int trials = 0;
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, 
        prune, 
        radius, 
        regrafts);
    for (auto regraft: regrafts) {
      trials++;
      if (testPruning(prune, regraft)) {
        Logger::timed << "\tFOUND VERY BETTER (LL=" 
            << _bestRecLL << ", hash=" << 
            _speciesTree->getHash() << " wrong_clades=" 
            << _unsupportedCladesNumber() << ")"<< std::endl;
        Logger::info << _speciesTree->getNode(prune)->label << " " << _speciesTree->getNode(regraft)->label << std::endl;
        return veryLocalSearch(prune);
      }
    }
  }
  return _bestRecLL;
}

double SpeciesTreeOptimizer::transferSearch()
{
  auto bestLL = computeRecLikelihood();
  Logger::timed << "[Species tree search, transfer-guided]" << " Starting species transfer search, bestLL=" 
    << bestLL <<  ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() 
    << ")" <<std::endl;
  double newLL = bestLL;
  MovesBlackList blacklist;
  unsigned int index = 0;
  bool maxImprovementsReached = true;
  do {
    bestLL = newLL;
    if (maxImprovementsReached) {
      bestLL = _bestRecLL = optimizeDTLRates();
    }
    newLL = transferRound(blacklist, maxImprovementsReached);
    Logger::info << "  Accepted: " << _okForClades << std::endl;
    Logger::info << "  Rejected: " << _koForClades << std::endl;
    index++;
  } while (newLL - bestLL > 1.0);
  if (index == 1) {
      newLL = bestLL = _bestRecLL = optimizeDTLRates();
  }
  Logger::timed << "After transfer search: " << newLL << std::endl;
  saveCurrentSpeciesTreeId();
  return newLL;
}

double SpeciesTreeOptimizer::sprSearch(unsigned int radius)
{
  double bestLL = computeRecLikelihood();
  Logger::timed << "[Species tree search, SPR moves]" << " Starting species SPR search, radius=" 
    << radius << ", bestLL=" << bestLL  << ", hash=" << _speciesTree->getHash() << " wrong_clades=" << _unsupportedCladesNumber() << ")" <<std::endl;
  double newLL = bestLL;
  do {
    bestLL = newLL;
    newLL = fastSPRRound(radius);
  } while (newLL - bestLL > 0.001);
  Logger::timed << "After normal search: " << bestLL << std::endl;
  saveCurrentSpeciesTreeId();
  Logger::info << "  Accepted: " << _okForClades << std::endl;
  Logger::info << "  Rejected: " << _koForClades << std::endl;
  return newLL;
}
  
ModelParameters SpeciesTreeOptimizer::computeOptimizedRates() 
{
  if (_userDTLRates || _optimizationCriteria == SupportedClades) {
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
  if (_userDTLRates || _optimizationCriteria == SupportedClades) {
    return computeRecLikelihood();
  }
  Logger::timed << "optimize rates " << std::endl;
  _modelRates = computeOptimizedRates();
  unsigned int i = 0;
  for (auto &evaluation: _evaluations) {
    evaluation->setRates(_modelRates.getRates(i++));
  }
  Logger::timed << "optimize rates done (LL=" << computeRecLikelihood() << ")" << std::endl;
  if (!_modelRates.info.perFamilyRates) {
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


double SpeciesTreeOptimizer::computeRecLikelihood()
{
  double res = 0.0;
  switch (_optimizationCriteria) {
  case ReconciliationLikelihood: 
    for (auto &evaluation: _evaluations) {
      res += evaluation->evaluate();
    }
    ParallelContext::sumDouble(res);
    break;
  case SupportedClades:
    res = -static_cast<double>(_unsupportedCladesNumber());
    break;
  }
  return res;

}

void SpeciesTreeOptimizer::newBestTreeCallback()
{
  saveCurrentSpeciesTreeId();
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
    _evaluations[i] = std::make_shared<ReconciliationEvaluation>(_speciesTree->getTree(), *tree.geneTree, tree.mapping, _modelRates.info);
    _evaluations[i]->setRates(_modelRates.getRates(i));
    _evaluations[i]->setPartialLikelihoodMode(PartialLikelihoodMode::PartialSpecies);
  }
  _previousGeneRoots.resize(_evaluations.size());
  std::fill(_previousGeneRoots.begin(), _previousGeneRoots.end(), nullptr);
}
  
void SpeciesTreeOptimizer::beforeTestCallback()
{
  if (_modelRates.info.rootedGeneTree) {
    for (unsigned int i = 0; i < _evaluations.size(); ++i) {
      _previousGeneRoots[i] = _evaluations[i]->getRoot();
    }
  }
}

void SpeciesTreeOptimizer::rollbackCallback()
{
  if (_modelRates.info.rootedGeneTree) {
    for (unsigned int i = 0; i < _evaluations.size(); ++i) {
      _evaluations[i]->setRoot(_previousGeneRoots[i]);
    }
  }
}
  
void SpeciesTreeOptimizer::optimizeGeneRoots()
{
  for (unsigned int i = 0; i < _evaluations.size(); ++i) {
    _evaluations[i]->setRoot(nullptr);
  }
  
}

void SpeciesTreeOptimizer::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  for (auto &evaluation: _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
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
  ParallelContext::barrier();
  std::remove(getCladesSetPath(_outputDir, ParallelContext::getRank()).c_str());
  Logger::timed << "Number number of supported bipartitions: " << _geneClades.size()  << std::endl;
}


unsigned int SpeciesTreeOptimizer::_unsupportedCladesNumber()
{
  auto speciesClades = Clade::buildCladeSet(_speciesTree->getTree());
  unsigned int intersectionSize =  intersection_size(speciesClades, 
      _geneClades);
  return speciesClades.size() - intersectionSize;
}
  
