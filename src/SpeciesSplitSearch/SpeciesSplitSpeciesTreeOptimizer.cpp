#include "SpeciesSplitSpeciesTreeOptimizer.hpp" 
#include <ccp/RootedSpeciesSplitScore.hpp>
#include <ccp/UnrootedSpeciesSplitScore.hpp>
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <util/Paths.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}
  
SpeciesSplitEvaluator::SpeciesSplitEvaluator(PLLRootedTree &speciesTree,
    const Families &families,
    bool rootedScore):
    _speciesTree(&speciesTree),
    _families(&families),
    _splits(speciesTree.getLabels(true), rootedScore)
{
  for (auto &family: families) {
    PLLUnrootedTree geneTree(family.startingGeneTree);
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    _splits.addGeneTree(geneTree, mapping);
  }
  _splits.computeVector();
  if (rootedScore) {
    _score = std::make_unique<RootedSpeciesSplitScore>(*_speciesTree, _splits);
  } else {
    _score = std::make_unique<UnrootedSpeciesSplitScore>(*_speciesTree, _splits);
  }
}

double SpeciesSplitEvaluator::computeLikelihood()
{
  return computeLikelihoodFast();
}

double SpeciesSplitEvaluator::computeLikelihoodFast()
{
  _score->updateSpeciesTree(*_speciesTree);
  auto res = _score->getScore();
  return res;
}


  
void SpeciesSplitEvaluator::getTransferInformation(PLLRootedTree &speciesTree,
    TransferFrequencies &transferFrequencies,
    PerSpeciesEvents &perSpeciesEvents)
{
  // this is duplicated code from Routines...
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int labelsNumber = idToLabel.size();
  const VectorUint zeros(labelsNumber, 0);
  transferFrequencies.count = MatrixUint(labelsNumber, zeros);
  transferFrequencies.idToLabel = idToLabel;
  perSpeciesEvents = PerSpeciesEvents(speciesTree.getNodesNumber());
  RecModelInfo info;
  info.pruneSpeciesTree = false;
  PerCoreGeneTrees geneTrees(*_families);
  for (const auto &geneTree: geneTrees.getTrees()) {
    auto &family = (*_families)[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    ReconciliationEvaluation evaluation(
        speciesTree,
        *geneTree.geneTree,
        mapping,
        info);
    Scenario scenario;
    Parameters rates(0.1, 0.1, 0.1);
    evaluation.setRates(rates);
    evaluation.evaluate();
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

void SpeciesSplitEvaluator::fillPerFamilyLikelihoods(
    PerFamLL &perFamLL)
{
}


SpeciesSplitSpeciesTreeOptimizer::SpeciesSplitSpeciesTreeOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    const std::string &outputDir,
    bool rootedScore):
  _speciesTree(std::make_unique<SpeciesTree>(speciesTreeFile)),
  _evaluator(_speciesTree->getTree(), families, rootedScore),
  _outputDir(outputDir),
  _searchState(*_speciesTree,
      Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick")),
  _rootedScore(rootedScore)
{
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
  _speciesTree->addListener(this);
  ParallelContext::barrier();
  Logger::timed << "Initial ll=" << _evaluator.computeLikelihood() << std::endl;
}

void SpeciesSplitSpeciesTreeOptimizer::optimize()
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
  do {
    rootSearch(3);
    if (index++ % 2 == 0) {
      transferSearch();
    } else {
      sprSearch(3);
    }
    hash1 = _speciesTree->getHash();
  }
  while(testAndSwap(hash1, hash2));
}

double SpeciesSplitSpeciesTreeOptimizer::rootSearch(unsigned int maxDepth)
{
  if (!_rootedScore) {
    return _searchState.bestLL;
  }
  Logger::info << std::endl;
  Logger::timed << "[Species search] Root search with depth=" << maxDepth << std::endl;
  RootLikelihoods rootLikelihoods;
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      _evaluator,
      _searchState,
      maxDepth,
      &rootLikelihoods,
      nullptr);
  saveCurrentSpeciesTreeId();
  return _searchState.bestLL;
}

double SpeciesSplitSpeciesTreeOptimizer::sprSearch(unsigned int radius)
{
  SpeciesSPRSearch::SPRSearch(*_speciesTree,
      _evaluator,
      _searchState,
      radius);
  Logger::timed << "After normal search: LL=" << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}


void SpeciesSplitSpeciesTreeOptimizer::onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate)
{
}



std::string SpeciesSplitSpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  return res;
}

void SpeciesSplitSpeciesTreeOptimizer::saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly)
{
  _speciesTree->saveToFile(str, masterRankOnly);
  if (masterRankOnly) {
    ParallelContext::barrier();
  }
}
  
double SpeciesSplitSpeciesTreeOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
    *_speciesTree,
    _evaluator,
    _searchState);
  Logger::timed << "After normal search: LL=" 
    << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}
  


