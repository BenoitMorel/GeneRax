#include "MiniBMEOptimizer.hpp" 
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <util/Paths.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <corax/corax.h>

double USearchMiniBMEEvaluator::eval(PLLUnrootedTree &tree)
{
  _lastScore = -_miniBME->computeBME(tree);
  //corax_unode_t *fake;
  //_miniBME.getBestSPR(tree, fake, fake);
  return _lastScore;
}

  
bool USearchMiniBMEEvaluator::computeAndApplyBestSPR(PLLUnrootedTree &tree,
    unsigned int maxRadiusWithoutImprovement)
{
  _lastScore = -_miniBME->computeBME(tree);
  double diff = 0.0;
  corax_unode_t *bestPruneNode = nullptr;
  corax_unode_t *bestRegraftNode = nullptr;
  _miniBME->getBestSPR(tree, 
      maxRadiusWithoutImprovement,
      bestPruneNode, 
      bestRegraftNode,
      diff);
  double epsilon = 0.00000001;
  if (diff <= epsilon) {
    Logger::info << "Local search failed, trying with max radius..." << std::endl;
    _miniBME->getBestSPR(tree, 
        99999999,
        bestPruneNode, 
        bestRegraftNode,
        diff);
  }
  if (diff > epsilon) {
    assert(bestPruneNode);
    assert(bestRegraftNode);
    auto ok = corax_utree_spr(bestPruneNode->back, 
        bestRegraftNode, 
        nullptr);
    assert(ok);
    Logger::timed << "Score: " << _lastScore + diff << std::endl;
    /*
    Logger::timed << "Estimated diff " << diff << std::endl; 
    Logger::timed << "Real      diff " << -_miniBME.computeBME(tree) - _lastScore << std::endl; 
    */
    return true;
  } else {
    return false;
  }
}


MiniBMEOptimizer::MiniBMEOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    double minbl,
    bool missingData,
    const std::string &outputDir):
  _speciesTree(std::make_unique<SpeciesTree>(speciesTreeFile)),
  _outputDir(outputDir),
  _minbl(minbl),
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
    _minbl,
    _missingData);
  Logger::info << "Starting score: " << evaluator.eval(speciesTree) << std::endl;
  bool ok = true;
  unsigned int it = 0;
  while (ok) {
    ok = evaluator.computeAndApplyBestSPR(speciesTree, 5);
    ++it;
  }
  Logger::info << "SPR search stopped after " << it << " iterations" << std::endl;
  Logger::info << "Final score: " << evaluator.eval(speciesTree) << std::endl;
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
  


