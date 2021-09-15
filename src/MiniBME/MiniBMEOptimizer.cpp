#include "MiniBMEOptimizer.hpp" 
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <util/Paths.hpp>
#include <search/SpeciesSPRSearch.hpp>
#include <search/SpeciesTransferSearch.hpp>
#include <search/UNNISearch.hpp>
#include <likelihoods/reconciliation_models/UndatedDTLModel.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
  

double USearchMiniBMEEvaluator::eval(PLLUnrootedTree &tree)
{
  _lastScore = -_miniBME.computeBME(tree);
  return _lastScore;
}

double USearchMiniBMEEvaluator::evalNNI(PLLUnrootedTree &tree,
    UNNIMove &move)
{
  /*  
  auto before = eval(tree);
  auto diff2 = _miniBME.computeNNIDiff(tree, move);
  move.apply();
  auto after = eval(tree);
  move.apply(); // rollback
  auto diff1 = before - after;
  std::cerr << "diffs: " << diff1 << " " << diff2 << " " << diff1/diff2 << std::endl;
  return after;
  */
  return _lastScore - _miniBME.computeNNIDiff(tree, move);
}

static bool testAndSwap(size_t &hash1, size_t &hash2) {
  std::swap(hash1, hash2);
  return hash1 != hash2;
}
  
MiniBMEEvaluator::MiniBMEEvaluator(PLLRootedTree &speciesTree,
    const Families &families,
    bool missingData):
    _speciesTree(&speciesTree),
    _families(&families),
    _miniBME(speciesTree, families, missingData)
{
}

double MiniBMEEvaluator::computeLikelihood()
{
  return computeLikelihoodFast();
}

double MiniBMEEvaluator::computeLikelihoodFast()
{
  PLLUnrootedTree speciesTree(*_speciesTree);
  auto res = -_miniBME.computeBME(speciesTree);
  return res;
}


  
void MiniBMEEvaluator::getTransferInformation(PLLRootedTree &speciesTree,
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

void MiniBMEEvaluator::fillPerFamilyLikelihoods(
    PerFamLL &perFamLL)
{
}


MiniBMEOptimizer::MiniBMEOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    bool missingData,
    const std::string &outputDir):
  _speciesTree(std::make_unique<SpeciesTree>(speciesTreeFile)),
  _evaluator(_speciesTree->getTree(), families, missingData),
  _outputDir(outputDir),
  _missingData(missingData),
  _families(families),
  _searchState(*_speciesTree,
      Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"))
{
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
  _speciesTree->addListener(this);
  ParallelContext::barrier();
  Logger::timed << "Initial ll=" << _evaluator.computeLikelihood() << std::endl;
}

void MiniBMEOptimizer::optimize()
{
  if (_missingData) {
    _searchState.bestLL = _evaluator.computeLikelihood();
    //sprSearch(3);
    size_t hash1 = 0;
    size_t hash2 = 0;
    unsigned int index = 0;
    do {
      if (index++ % 2 == 0) {
        transferSearch();
      } else {
        sprSearch(3);
      }
      hash1 = _speciesTree->getHash();
    }
    while(testAndSwap(hash1, hash2));
  } else {
    PLLUnrootedTree speciesTree(_speciesTree->getTree());
    USearchMiniBMEEvaluator evaluator(speciesTree,
      _families,
      _missingData);
    UNNISearch search(speciesTree, evaluator);
    search.search();
    speciesTree.save(Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"));
  }
}


double MiniBMEOptimizer::sprSearch(unsigned int radius)
{
  SpeciesSPRSearch::SPRSearch(*_speciesTree,
      _evaluator,
      _searchState,
      radius);
  Logger::timed << "After normal search: LL=" << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}


void MiniBMEOptimizer::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
}



std::string MiniBMEOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  return res;
}

void MiniBMEOptimizer::saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly)
{
  _speciesTree->saveToFile(str, masterRankOnly);
  if (masterRankOnly) {
    ParallelContext::barrier();
  }
}
  
double MiniBMEOptimizer::transferSearch()
{
  SpeciesTransferSearch::transferSearch(
    *_speciesTree,
    _evaluator,
    _searchState);
  Logger::timed << "After normal search: LL=" 
    << _searchState.bestLL << std::endl;
  return _searchState.bestLL;
}
  


