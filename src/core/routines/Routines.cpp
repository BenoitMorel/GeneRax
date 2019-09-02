
#include "Routines.hpp"

#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/FileSystem.hpp>
#include <likelihoods/LibpllEvaluation.hpp>

void Routines::optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    RecModel recModel,
    Families &families,
    bool perSpeciesRates, 
    DTLRatesVector &rates,
    long &sumElapsed) 
{
  if (userDTLRates) {
    return;
  }
  auto start = Logger::getElapsedSec();
  Logger::timed << "Start optimizing rates..." << std::endl;
  PerCoreGeneTrees geneTrees(families);
  bool ok = geneTrees.checkMappings(speciesTreeFile); 
  if (!ok) {
    Logger::info << "INVALID MAPPINGS" << std::endl;
    ParallelContext::abort(42);
  }
  pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
  if (!perSpeciesRates) {
    rates = DTLOptimizer::optimizeDTLRates(geneTrees, speciesTree, recModel);
  } else {
    rates = DTLOptimizer::optimizeDTLRatesVector(geneTrees, speciesTree, recModel);
  }
  pll_rtree_destroy(speciesTree, 0);
  ParallelContext::barrier(); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "Finished optimizing rates: "
    << "Loglk=" << rates.getLL() 
    << " (after " << elapsed << "s)" << std::endl;
}

void Routines::inferReconciliation(
    const std::string &speciesTreeFile,
    Families &families,
    RecModel model,
    DTLRatesVector &rates,
    const std::string &outputDir
    )
{
  pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
  PerCoreGeneTrees geneTrees(families);
  std::string reconciliationsDir = FileSystem::joinPaths(outputDir, "reconciliations");
  FileSystem::mkdir(reconciliationsDir, true);
  auto speciesNodesCount = speciesTree->tip_count + speciesTree->inner_count;
  std::vector<double> dup_count(speciesNodesCount, 0.0);
  ParallelContext::barrier();
  for (auto &tree: geneTrees.getTrees()) {
    std::string eventCountsFile = FileSystem::joinPaths(reconciliationsDir, tree.name + "_eventCounts.txt");
    std::string treeWithEventsFileNHX = FileSystem::joinPaths(reconciliationsDir, tree.name + "_reconciliated.nhx");
    std::string treeWithEventsFileRecPhyloXML = FileSystem::joinPaths(reconciliationsDir, tree.name + "_reconciliated.recml");
    Scenario scenario;
    ReconciliationEvaluation evaluation(speciesTree, tree.mapping, model, true);
    evaluation.setRates(rates);
    evaluation.evaluate(tree.tree);
    evaluation.inferMLScenario(scenario);
    scenario.saveEventsCounts(eventCountsFile, false);
    scenario.saveReconciliations(treeWithEventsFileNHX, NHX, false);
    scenario.saveReconciliations(treeWithEventsFileRecPhyloXML, RecPhyloXML, false);
    for (auto &event: scenario.getEvents()) {
      if (event.type == EVENT_D) {
        dup_count[event.speciesNode]++;
      } 
    }
  }
  ParallelContext::sumVectorDouble(dup_count);
  auto initialSpeciesTreeStr = pll_rtree_export_newick(speciesTree->root, 0);
  free(initialSpeciesTreeStr);
  auto speciesTreeStr = pll_rtree_export_newick(speciesTree->root, 0);
  free(speciesTreeStr);
  pll_rtree_destroy(speciesTree, 0);
}

bool Routines::createRandomTrees(const std::string &geneRaxOutputDir, Families &families)
{
  std::string startingTreesDir = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
  bool startingTreesDirCreated = false;
  for (auto &family: families) {
    if (family.startingGeneTree == "__random__") {
        if (!startingTreesDirCreated) {
          FileSystem::mkdir(startingTreesDir, true);
          startingTreesDirCreated = true;
        } 
        family.startingGeneTree = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
        family.startingGeneTree = FileSystem::joinPaths(family.startingGeneTree, family.name + ".newick");
        if (ParallelContext::getRank() == 0) {
          LibpllEvaluation::createAndSaveRandomTree(family.alignmentFile, family.libpllModel, family.startingGeneTree);
        }
    }
  }
  ParallelContext::barrier();
  return startingTreesDirCreated;
}

void Routines::gatherLikelihoods(Families &families,
    double &totalLibpllLL,
    double &totalRecLL)
{
  totalRecLL = 0.0;
  totalLibpllLL = 0.0;
  unsigned int familiesNumber = static_cast<unsigned int>(families.size());
  for (auto i = ParallelContext::getBegin(familiesNumber); i < ParallelContext::getEnd(familiesNumber); ++i) {
    auto &family = families[i];
    std::ifstream is(family.statsFile);
    double libpllLL = 0.0;
    double recLL = 0.0;
    is >> libpllLL;
    is >> recLL;
    totalRecLL += recLL;
    totalLibpllLL += libpllLL;
  }
  ParallelContext::sumDouble(totalRecLL);
  ParallelContext::sumDouble(totalLibpllLL);
}

