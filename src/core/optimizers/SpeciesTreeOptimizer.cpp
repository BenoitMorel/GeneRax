#include "SpeciesTreeOptimizer.hpp"

#include <optimizers/DTLOptimizer.hpp>
#include <IO/FileSystem.hpp>
#include <routines/GeneRaxMaster.hpp>
#include <routines/Routines.hpp>
#include <algorithm>

#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>

static std::string getStepTag(bool fastMove)
{
  static std::string fastMoveString("[Species tree search - Fast moves]");
  static std::string slowMoveString("[Species tree search - Slow moves]");
  return fastMove ? fastMoveString : slowMoveString;
}

SpeciesTreeOptimizer::SpeciesTreeOptimizer(const std::string speciesTreeFile, 
    const Families &initialFamilies, 
    RecModel model,
    double supportThreshold,
    const std::string &outputDir,
    const std::string &execPath):
  _speciesTree(nullptr),
  _geneTrees(nullptr),
  _initialFamilies(initialFamilies),
  _currentFamilies(initialFamilies),
  _recModel(model),
  _outputDir(outputDir),
  _execPath(execPath),
  _geneTreeIteration(1000000000), // we need to find a better way for avoiding directories collision
  _supportThreshold(supportThreshold),
  _lastRecLL(-std::numeric_limits<double>::infinity()),
  _lastLibpllLL(-std::numeric_limits<double>::infinity()),
  _bestRecLL(-std::numeric_limits<double>::infinity()),
  _bestLibpllLL(-std::numeric_limits<double>::infinity())
  
{
  if (speciesTreeFile == "random") {
    _speciesTree = std::make_unique<SpeciesTree>(initialFamilies);
    _speciesTree->setGlobalRates(Parameters(0.1, 0.2, 0.1));
    setGeneTreesFromFamilies(initialFamilies);
  } else {
    _speciesTree = std::make_unique<SpeciesTree>(speciesTreeFile);
    _speciesTree->setGlobalRates(Parameters(0.1, 0.2, 0.1));
    setGeneTreesFromFamilies(initialFamilies);
    optimizeDTLRates();
  }
  _speciesTree->saveToFile(FileSystem::joinPaths(_outputDir, "starting_species_tree.newick"), true);
  _speciesTree->addListener(this);
  std::string subsamplesPath = FileSystem::joinPaths(_outputDir, "subsamples");
  FileSystem::mkdir(FileSystem::joinPaths(_outputDir, "sub_genes_opt"), true);
  FileSystem::mkdir(subsamplesPath, true);
  saveCurrentSpeciesTreeId();
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

void SpeciesTreeOptimizer::rootExhaustiveSearch(bool doOptimizeGeneTrees)
{
  //Logger::timed << getStepTag(!doOptimizeGeneTrees) << " Trying to re-root the species tree" << std::endl;

  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  unsigned int geneRadius = doOptimizeGeneTrees ? 1 : 0;
  double bestLL = computeLikelihood(geneRadius);
  unsigned int visits = 1;
  movesHistory.push_back(0);
  rootExhaustiveSearchAux(*_speciesTree, 
      *_geneTrees, 
      _recModel, 
      doOptimizeGeneTrees, 
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits); 
  movesHistory[0] = 1;
  rootExhaustiveSearchAux(*_speciesTree, 
      *_geneTrees, 
      _recModel, 
      doOptimizeGeneTrees, 
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits); 
  assert (visits == 2 * _speciesTree->getTree().getLeavesNumber() - 3);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(*_speciesTree, bestMovesHistory[i]);
  }
}
 



double SpeciesTreeOptimizer::fastSPRRound(unsigned int radius)
{
  bool check = false;
  bool tryApproxFirst = Enums::implementsApproxLikelihood(_recModel);
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, prunes);
  _bestRecLL = computeRecLikelihood();
  auto refApproxLL = computeApproxRecLikelihood();
  assert (fabs(_bestRecLL - refApproxLL) < 0.01);
  auto hash1 = _speciesTree->getNodeIndexHash(); 
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      // Apply the move
      auto rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
      bool canTestMove = true;
      // Discard bad moves with an approximation of the likelihood function
      double approxRecLL;
      bool needFullRollback = false;
      if (tryApproxFirst) {
        approxRecLL = computeApproxRecLikelihood();
        //Logger::info << approxRecLL << std::endl;
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
          Logger::info << "new best tree " << _bestRecLL << " -> " << _lastRecLL << std::endl;
          newBestTreeCallback();
          return _lastRecLL;
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
      saveCurrentSpeciesTreeId();
      newBestTreeCallback();
      return newBestLL;
    }
    SpeciesTreeOperator::reverseSPRMove(*_speciesTree, em.prune, rollback);
  } 
  return bestLL;
}

double SpeciesTreeOptimizer::sprSearch(unsigned int radius, bool doOptimizeGeneTrees)
{
  unsigned int geneRadius = doOptimizeGeneTrees ? 1 : 0;
  double bestLL = doOptimizeGeneTrees ? computeLikelihood(geneRadius) : computeRecLikelihood();
  Logger::timed << getStepTag(!doOptimizeGeneTrees) << " Starting species SPR search, radius=" 
    << radius << ", bestLL=" << bestLL << ")" <<std::endl;
  double newLL = bestLL;
  do {
    bestLL = newLL;
    if (doOptimizeGeneTrees) {
      newLL = slowSPRRound(radius, bestLL); 
    } else {
      newLL = fastSPRRound(radius);
    }
  } while (newLL - bestLL > 0.001);
  saveCurrentSpeciesTreeId();
  return newLL;
}
  
Parameters SpeciesTreeOptimizer::computeOptimizedRates() const
{
  return DTLOptimizer::optimizeParametersGlobalDTL(*_geneTrees, _speciesTree->getTree(), _recModel);
}
  
void SpeciesTreeOptimizer::optimizeDTLRates()
{
  _speciesTree->setGlobalRates(computeOptimizedRates());
  for (auto &evaluation: _evaluations) {
    evaluation->setRates(_speciesTree->getRatesVector());
  }
}
  
void SpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  saveCurrentSpeciesTreePath(FileSystem::joinPaths(_outputDir, name), masterRankOnly);
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
  auto rates = _speciesTree->getRatesVector();
  std::string resultName = "proposals";
  unsigned int iterationsNumber = 1;
  bool inPlace = false; 
  bool perFamilyDTLRates = false;
  if (radius == 1) {
    iterationsNumber = 2;
  }
  for (unsigned i = 0; i < iterationsNumber; ++i) {
    Logger::mute();
    GeneRaxMaster::optimizeGeneTrees(_currentFamilies, 
      _recModel, rates, _outputDir, resultName, 
      _execPath, speciesTree, recOpt, perFamilyDTLRates, rootedGeneTree, 
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
  return ll;
}

double SpeciesTreeOptimizer::computeApproxRecLikelihood()
{
  double ll = 0.0;
  for (auto &evaluation: _evaluations) {
    ll += evaluation->evaluate(true);
  }
  ParallelContext::sumDouble(ll);
  return ll;
}

void SpeciesTreeOptimizer::newBestTreeCallback()
{
  _bestLibpllLL = _lastLibpllLL;
  _bestRecLL = _lastRecLL;
}
  
void SpeciesTreeOptimizer::setGeneTreesFromFamilies(const Families &families)
{
  _geneTrees = std::make_unique<PerCoreGeneTrees>(families);
  updateEvaluations();
}
  
void SpeciesTreeOptimizer::updateEvaluations()
{
  assert(_geneTrees);
  auto &trees = _geneTrees->getTrees();
  _evaluations.resize(trees.size());
  for (unsigned int i = 0; i < trees.size(); ++i) {
    auto &tree = trees[i];
    _evaluations[i] = std::make_shared<ReconciliationEvaluation>(_speciesTree->getTree(), *tree.geneTree, tree.mapping, _recModel, false);
    _evaluations[i]->setRates(_speciesTree->getRatesVector());
    _evaluations[i]->setPartialLikelihoodMode(PartialLikelihoodMode::PartialSpecies);
  }
}

void SpeciesTreeOptimizer::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  for (auto &evaluation: _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}



