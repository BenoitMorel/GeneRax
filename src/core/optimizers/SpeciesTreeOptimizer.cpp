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
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>

static std::unique_ptr<SpeciesTree> makeSpeciesTree(const std::string &speciesTreeFile,
    const Families &initialFamilies)
{
  if (speciesTreeFile == "random") {
    return std::make_unique<SpeciesTree>(initialFamilies);
  } else {
    return std::make_unique<SpeciesTree>(speciesTreeFile);
  }

}

SpeciesTreeOptimizer::SpeciesTreeOptimizer(const std::string speciesTreeFile, 
    const Families &initialFamilies, 
    const RecModelInfo &recModelInfo,
    const Parameters &startingRates,
    bool userDTLRates,
    const std::string &outputDir,
    const SpeciesTreeSearchParams &searchParams):
  _speciesTree(makeSpeciesTree(speciesTreeFile,
        initialFamilies)),
  _geneTrees(nullptr),
  _initialFamilies(initialFamilies),
  _outputDir(outputDir),
  _firstOptimizeRatesCall(true),
  _userDTLRates(userDTLRates),
  _modelRates(startingRates, 1, recModelInfo),
  _searchParams(searchParams),
  _okForClades(0),
  _koForClades(0),
  _searchState(*_speciesTree, 
      Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"))
{

  _modelRates.info.perFamilyRates = false; // we set it back a few
                                           // lines later
  setGeneTreesFromFamilies(initialFamilies);
  _modelRates = ModelParameters(startingRates, 
      _geneTrees->getTrees().size(),
      recModelInfo);
  _speciesTree->addListener(this);
  saveCurrentSpeciesTreeId();
  _computeAllGeneClades();
  _searchState.farFromPlausible &= 
    _unsupportedCladesNumber() > std::max<unsigned int>(
      _speciesTree->getTree().getNodesNumber() / 4, 1);
  if (_searchState.farFromPlausible) {
    Logger::info << "The starting tree is far from being plausible" << std::endl;
  }
}

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}

void SpeciesTreeOptimizer::optimize(SpeciesSearchStrategy strategy)
{
  computeRecLikelihood();
  size_t hash1 = 0;
  size_t hash2 = 0;
  unsigned int index = 0;
  switch (strategy) {
  case SpeciesSearchStrategy::SPR:
    for (unsigned int radius = 1; radius <= _searchParams.sprRadius; ++radius) {
      optimizeDTLRates();
      sprSearch(radius);
    }
    break;
  case SpeciesSearchStrategy::TRANSFERS:
    transferSearch();
    rootSearch(_searchParams.rootBigRadius, false);
    transferSearch();
    rootSearch(_searchParams.rootBigRadius, true);
    break;
  case SpeciesSearchStrategy::HYBRID:
    /**
     *  Alternate transfer search and normal
     *  SPR search, until one does not find
     *  a better tree. Run each at least once.
     */
    if (!_searchState.farFromPlausible) {
      optimizeDTLRates();
      computeRecLikelihood();
      rootSearch(_searchParams.rootSmallRadius, false);
    }
    do {
      if (index++ % 2 == 0) {
        transferSearch();
      } else {
        sprSearch(_searchParams.sprRadius);
      }
      if (!_searchState.farFromPlausible) {
        rootSearch(_searchParams.rootSmallRadius, false);
      }
      hash1 = _speciesTree->getHash();
    }
    while(testAndSwap(hash1, hash2));
    rootSearch(_searchParams.rootBigRadius, true);
    break;
  case SpeciesSearchStrategy::REROOT:
    rootSearch(_searchParams.rootBigRadius, true);
    break;
  case SpeciesSearchStrategy::EVAL:
    optimizeDTLRates(true);
    Logger::info << "Reconciliation likelihood: " << computeRecLikelihood() << std::endl;
    break;
  case SpeciesSearchStrategy::SKIP:
    assert(false);
  }
}


SpeciesTreeOptimizer::~SpeciesTreeOptimizer()
{
  _speciesTree->removeListener(this);
}
  
double SpeciesTreeOptimizer::rootSearch(unsigned int maxDepth,
    bool outputConsel)
{
  Logger::info << std::endl;
  Logger::timed << "[Species search] Root search with depth=" << maxDepth << std::endl;
  TreePerFamLLVec treePerFamLLVec;  
  RootLikelihoods rootLikelihoods;
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      _evaluator,
      _searchState,
      maxDepth,
      &rootLikelihoods,
      (outputConsel ? &treePerFamLLVec : nullptr)
      );
  saveCurrentSpeciesTreeId();
  {
    auto newick = _speciesTree->getTree().getNewickString();
    PLLRootedTree tree(newick, false); 
    rootLikelihoods.fillTree(tree);
    auto out = Paths::getSpeciesTreeFile(_outputDir, 
        "species_tree_llr.newick");
    tree.save(out);
  }
  if (outputConsel) {
    std::string treesOutput = Paths::getConselTreeList(_outputDir, 
        "roots"); 
    std::string llOutput  = Paths::getConselLikelihoods(_outputDir, 
        "roots"); 
    Logger::info << "Saving per-family likelihoods into: " 
      << llOutput << std::endl;
    Logger::info << "Saving the corresponding trees into: "
      << treesOutput << std::endl;
    savePerFamilyLikelihoods(treePerFamLLVec,
      treesOutput,
      llOutput);
    
  }
  auto ll = computeRecLikelihood();
  Logger::timed << "[Species search] After root search: LL=" << ll << std::endl;
  return ll;
}


std::vector<double> SpeciesTreeOptimizer::_getSupport()
{
  std::string temp = FileSystem::joinPaths(_outputDir, "tmp");
  std::vector<double> idToSupport;
  bool paralogyAware = true;
  unsigned int eqpicRadius = 3;
  ICCalculator::computeScores(_speciesTree->getTree(),
      _initialFamilies,
      paralogyAware,
      eqpicRadius,
      temp,
      idToSupport);
  for (auto node: _speciesTree->getTree().getLeaves()) {
    idToSupport[node->node_index] = 1.0;
  }
  return idToSupport;
}



double SpeciesTreeOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
      *_speciesTree,
    _evaluator,
    _searchState);
  Logger::timed << "After normal search: LL=" << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}

double SpeciesTreeOptimizer::sprSearch(unsigned int radius)
{
  SpeciesSPRSearch::SPRSearch(*_speciesTree,
    _evaluator,
    _searchState,
    radius);
  Logger::timed << "After normal search: LL=" 
    << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}
  
ModelParameters SpeciesTreeOptimizer::computeOptimizedRates(bool thorough) 
{
  if (_userDTLRates) {
    return _modelRates;
  }
  auto rates = _modelRates;
  OptimizationSettings settings;
  double ll = computeRecLikelihood();
  if (!thorough) {
    settings.lineSearchMinImprovement = 10.0;
    settings.minAlpha = 0.01;
    settings.optimizationMinImprovement = std::max(3.0, ll / 1000.0);
  }
  rates =  DTLOptimizer::optimizeModelParameters(_evaluations, !_firstOptimizeRatesCall, rates, settings);
  _firstOptimizeRatesCall = false;
  return rates;
}
  
double SpeciesTreeOptimizer::optimizeDTLRates(bool thorough)
{
  if (_userDTLRates) {
    return computeRecLikelihood();
  }
  //Logger::timed << "[Species search] Start rates optimization " << std::endl;
  _modelRates = computeOptimizedRates(thorough);
  unsigned int i = 0;
  for (auto &evaluation: _evaluations) {
    evaluation->setRates(_modelRates.getRates(i++));
  }
  //Logger::timed << "[Species search] Rates optimized! (LL=" << computeRecLikelihood() << ")" << std::endl;
  if (!_modelRates.info.perFamilyRates) {
    Logger::timed << "[Species search] Best rates: " << _modelRates.rates << std::endl;
  }
  return computeRecLikelihood();
}
  
std::string SpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
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
  for (auto &evaluation: _evaluations) {
    auto ll = evaluation->evaluate();
    res += ll;
  }
  ParallelContext::sumDouble(res);
  return res;

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
  _evaluator.init(_evaluations, 
      *_geneTrees,
      _modelRates,
      _modelRates.info.rootedGeneTree,
      _modelRates.info.pruneSpeciesTree);
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
}


unsigned int SpeciesTreeOptimizer::_unsupportedCladesNumber()
{
  auto speciesClades = Clade::buildCladeSet(_speciesTree->getTree());
  unsigned int intersectionSize =  intersection_size(speciesClades, 
      _geneClades);
  return speciesClades.size() - intersectionSize;
}

static std::string getSubtreeID(pll_rnode_t *subtree)
{
  if (!subtree->left) {
    return std::string(subtree->label);
  }
  std::string res("(");
  std::string id1 = getSubtreeID(subtree->left);
  std::string id2 = getSubtreeID(subtree->right);
  if  (id1 > id2) {
    std::swap(id1, id2);
  }
  return std::string("(") + id1 + "," + id2 + ")";
}
    
void SpeciesTreeOptimizer::savePerFamilyLikelihoods(
    const TreePerFamLLVec &treePerFamLLVec,
    const std::string &treesOutput,
    const std::string &llOutput)
{
  ParallelContext::barrier();
  if (ParallelContext::getRank() == 0 && treePerFamLLVec.size() != 0) {
    std::ofstream osLL(llOutput);
    std::ofstream osTrees(treesOutput);
    auto treesNumber = treePerFamLLVec.size();
    auto familiesNumber = treePerFamLLVec[0].second.size();
    osLL << treesNumber << " " << familiesNumber << std::endl;
    unsigned int index = 1;
    Logger::info << "trees = " << treePerFamLLVec.size() << std::endl;
    for (const auto &treePerFamLL: treePerFamLLVec) {
      const auto &tree = treePerFamLL.first;
      const auto perFamLL = treePerFamLL.second;
      std::string treeName = "tree";
      treeName += std::to_string(index++);
      osTrees << tree << std::endl;
      osLL << treeName;
      for (auto ll: perFamLL) {
        osLL << " " << ll;
      }
      osLL << std::endl;
    }
  }
  ParallelContext::barrier();

}


double SpeciesTreeLikelihoodEvaluator::computeLikelihood()
{ 
  if (_rootedGeneTrees) {
    for (auto evaluation: *_evaluations) {
      evaluation->setRoot(nullptr);
    }
  }
  return computeLikelihoodFast();
}

double SpeciesTreeLikelihoodEvaluator::computeLikelihoodFast()
{
  double sumLL = 0.0;
  for (auto &evaluation: *_evaluations) {
    auto ll = evaluation->evaluate();
    sumLL += ll;
  }
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}
  
bool SpeciesTreeLikelihoodEvaluator::providesFastLikelihoodImpl() const 
{
  return _rootedGeneTrees; 
}
  
void SpeciesTreeLikelihoodEvaluator::optimizeModelRates()
{
  Logger::info << "DTL RATES OPTIMIZATION IS NOT IMPLEMENTED" << std::endl;
}

void SpeciesTreeLikelihoodEvaluator::getTransferInformation(PLLRootedTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents)
{
  ParallelContext::barrier();
  unsigned int reconciliationSamples = 0; // use ML reconciliation
  Routines::getTransfersFrequencies(speciesTree,
    *_geneTrees,
    *_modelRates,
    reconciliationSamples,
    frequencies);
  const bool forceTransfers = true;
  Routines::getPerSpeciesEvents(speciesTree,
    *_geneTrees,
    *_modelRates,
    reconciliationSamples,
    perSpeciesEvents,
    forceTransfers);
}

void SpeciesTreeLikelihoodEvaluator::fillPerFamilyLikelihoods(
    PerFamLL &perFamLL)
{
  ParallelContext::barrier();
  std::vector<double> localLL;
  for (auto &evaluation: *_evaluations) {
    localLL.push_back(evaluation->evaluate());
  }
  ParallelContext::concatenateHetherogeneousDoubleVectors(
      localLL, perFamLL);
}

void SpeciesTreeLikelihoodEvaluator::pushRollback() 
{
  if (_rootedGeneTrees) {
    _previousGeneRoots.push(std::vector<pll_unode_t *>());
    for (auto evaluation: *_evaluations) {
      _previousGeneRoots.top().push_back(evaluation->getRoot());
    }
  }
}

void SpeciesTreeLikelihoodEvaluator::popAndApplyRollback() 
{
  if (_rootedGeneTrees) {
    for (unsigned int i = 0; i < _evaluations->size(); ++i) {
      (*_evaluations)[i]->setRoot(_previousGeneRoots.top()[i]);
    }
    _previousGeneRoots.pop();
  }
}
  



