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
  Logger::timed << "score=" << _lastScore << std::endl;
  return _lastScore;
}

static bool sprYeldsSameTree(corax_unode_t *p, corax_unode_t *r)
{
  assert(p);
  assert(r);
  assert(p->next);
  return (r == p) || (r == p->next) || (r == p->next->next)
    || (r == p->back) || (r == p->next->back) || (r == p->next->next->back);
}

bool isSPRMoveValid(PLLUnrootedTree &tree,
    corax_unode_t *prune, 
    corax_unode_t *regraft)
{
  // regraft should not be a child of prune
  auto pruneChildren = tree.getPostOrderNodesFrom(prune->back);
  for (auto child: pruneChildren) {
    if (regraft == child) {
      return false;
    }
    if (regraft->back == child) {
      return false;
    }
  }
  return !sprYeldsSameTree(prune, regraft);
}

bool isInScores(const std::vector<double> &scores, double score)
{
  for (auto s: scores) {
    if (fabs(s-score) < 0.000001) {
      return true;
    }
  }
  return false;
}

void addInvolvedNode(corax_unode_t *node, 
    std::unordered_set<corax_unode_t *> &involved)
{
  involved.insert(node);
  if (node->next) {
    involved.insert(node->next);
    involved.insert(node->next->next);
  }
}

bool wasInvolved(corax_unode_t *node,
    std::unordered_set<corax_unode_t *> &involved)
{
  return involved.find(node) != involved.end();
}
  
bool USearchMiniBMEEvaluator::computeAndApplyBestSPR(PLLUnrootedTree &tree,
    unsigned int maxRadiusWithoutImprovement)
{
  Logger::timed << "last score " << _lastScore << std::endl;
  std::vector<SPRMove> bestMoves;
  double epsilon = 0.00000001;
  _miniBME->getBestSPR(tree, 
      maxRadiusWithoutImprovement,
      bestMoves);
  if (bestMoves.size() == 0 || bestMoves[0].score < epsilon) {
    Logger::info << "Local search failed, trying with max radius..." << std::endl;
    _miniBME->getBestSPR(tree, 
        99999999,
        bestMoves);
  }
  bool better = false;
  unsigned int appliedMoves = 0;
  unsigned int maxAppliedMoves = 99999; // no max
  corax_tree_rollback_t emptyRollback;
  std::vector<corax_tree_rollback_t> rollbacks;
  std::vector<double> hackScore;
  std::unordered_set<corax_unode_t *> involved;
  double expectedDiff = 0.0;
  for (auto move: bestMoves) {
    if (isSPRMoveValid(tree, move.pruneNode->back, move.regraftNode) 
        && !isInScores(hackScore, move.score)
        && !wasInvolved(move.pruneNode, involved)
        && !wasInvolved(move.regraftNode, involved)
      ) {
      addInvolvedNode(move.pruneNode->back, involved);
      addInvolvedNode(move.regraftNode, involved);
      rollbacks.push_back(emptyRollback);   
      hackScore.push_back(move.score);
      /*
      Logger::info << move.score << " ";
      if (move.pruneNode->next) {
        Logger::info << "(" << move.pruneNode->label << "," << move.pruneNode->next->label << "," << move.pruneNode->next->next->label << ")";
      } else {
        Logger::info << move.pruneNode->label;
      }
      Logger::info  << " to ";
      if (move.regraftNode->next) {
        Logger::info << "(" << move.regraftNode->label << "," << move.regraftNode->next->label << "," << move.regraftNode->next->next->label << ")";
      } else {
        Logger::info << move.regraftNode->label;
      }
      Logger::info << std::endl;
      */
      expectedDiff += move.score;
      auto ok = corax_utree_spr(move.pruneNode->back, 
          move.regraftNode, 
          &rollbacks.back());
      appliedMoves++;
      assert(ok);
      better = true;
    }
    if (appliedMoves >= maxAppliedMoves) {
      break;
    }
    //Logger::timed << "Score: " << _lastScore + bestMoves[0].score << std::endl;
  }
  Logger::info << "Moves: " << appliedMoves << std::endl;
  double newScore = -_miniBME->computeBME(tree);
  double diff = newScore - _lastScore;
  if (appliedMoves) {
    Logger::info << "Expected: " << expectedDiff << " real:" << diff << " highest: " << bestMoves[0].score << std::endl;
  }
  if (newScore < _lastScore) {
    Logger::info << "New score " << newScore << " worse than last score " << _lastScore << std::endl;
    for (int i = rollbacks.size() - 1; i >= 1; --i) {
      corax_tree_rollback(&rollbacks[i]);
    }
    newScore = -_miniBME->computeBME(tree);
  }
  _lastScore = newScore;
    /*
    Logger::timed << "Estimated diff " << diff << std::endl; 
    Logger::timed << "Real      diff " << -_miniBME.computeBME(tree) - _lastScore << std::endl; 
    */
  return better;
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
  Logger::timed << "Evaluator initialization..." << std::endl;
  USearchMiniBMEEvaluator evaluator(speciesTree,
    _families,
    _minbl,
    _missingData);
  Logger::timed << "First score computation..." << std::endl;
  double startingScore = evaluator.eval(speciesTree); 
  Logger::info << "Starting score: " << startingScore << std::endl;
  bool ok = true;
  unsigned int it = 0;
  while (ok) {
    ok = evaluator.computeAndApplyBestSPR(speciesTree, 5);
    ++it;
  }
  Logger::info << "SPR search stopped after " << it << " iterations" << std::endl;
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
  


