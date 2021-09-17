#include "MiniBMEOptimizer.hpp" 
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <util/Paths.hpp>
#include <search/UNNISearch.hpp>
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
  return _lastScore - _miniBME.computeNNIDiff(tree, move);
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
  UNNISearch search(speciesTree, evaluator);
  search.search();
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
  


