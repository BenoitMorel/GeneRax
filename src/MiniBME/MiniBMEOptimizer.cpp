#include "MiniBMEOptimizer.hpp" 
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <util/Paths.hpp>
#include <search/UNNISearch.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <corax/corax.h>

double USearchMiniBMEEvaluator::eval(PLLUnrootedTree &tree)
{
  _lastScore = -_miniBME.computeBME(tree);
  //pll_unode_t *fake;
  //_miniBME.getBestSPR(tree, fake, fake);
  return _lastScore;
}

double USearchMiniBMEEvaluator::evalNNI(PLLUnrootedTree &tree,
    UNNIMove &move)
{
  /*
  auto before = eval(tree);
  Logger::info << tree.getUnrootedTreeHash() << " " << before << std::endl;
  auto diff2 = _miniBME.computeNNIDiff(tree, move);
  move.apply();
  auto after = eval(tree);
  Logger::info << tree.getUnrootedTreeHash() << " " << after << std::endl;
  move.apply(); // rollback
  auto diff1 = before - after;
  //Logger::info << "scores: " << before << " " << after << std::endl;
  //Logger::info << "diffs: " << diff1 << " " << diff2 << " " << diff1/diff2 << std::endl;
  eval(tree);
  return after;*/
  auto diff = _miniBME.computeNNIDiff(move);
  auto res = _lastScore - diff;
  Logger::info << "evalNNI diff=" << diff << std::endl;
  return res;
}
  
bool USearchMiniBMEEvaluator::computeAndApplyBestSPR(PLLUnrootedTree &tree)
{
  _lastScore = -_miniBME.computeBME(tree);
  double diff = 0.0;
  pll_unode_t *bestPruneNode = nullptr;
  pll_unode_t *bestRegraftNode = nullptr;
  _miniBME.getBestSPR(tree, 
      bestPruneNode, 
      bestRegraftNode,
      diff);
  if (diff > 0.00000001) {
    assert(bestPruneNode);
    assert(bestRegraftNode);
    auto ok = pllmod_utree_spr(bestPruneNode->back, 
        bestRegraftNode, 
        nullptr);
    assert(ok);
    Logger::timed << "Score: " << _lastScore + diff << std::endl;
    Logger::timed <<   "(diff=" << diff << ")" << std::endl;
    return true;
  } else {
    return false;
  }
}


MiniBMEOptimizer::MiniBMEOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    bool missingData,
    const std::string &outputDir):
  _speciesTree(std::make_unique<SpeciesTree>(speciesTreeFile)),
  _outputDir(outputDir),
  _missingData(missingData),
  _families(families)
{
  
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
  ParallelContext::barrier();
}

void MiniBMEOptimizer::optimize()
{
  PLLUnrootedTree speciesTree(_speciesTree->getTree());
  USearchMiniBMEEvaluator evaluator(speciesTree,
    _families,
    _missingData);
  /*
  UNNISearch search(speciesTree, evaluator);
  search.search();
  */
  bool ok = true;
  while (ok) {
    ok = evaluator.computeAndApplyBestSPR(speciesTree);
  }
  speciesTree.save(Paths::getSpeciesTreeFile(_outputDir, "inferred_species_tree.newick"));
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
  


