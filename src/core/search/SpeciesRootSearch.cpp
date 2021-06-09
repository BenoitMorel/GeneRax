
#include "SpeciesRootSearch.hpp"
#include <parallelization/ParallelContext.hpp>
#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>

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

void RootLikelihoods::saveValue(pll_rnode_t *subtree, double ll) 
{
  auto id = getSubtreeID(subtree); 
  idToLL[id] = ll;
}

void RootLikelihoods::fillTree(PLLRootedTree &tree)
{
  std::vector<double> nodeIdToLL(tree.getNodesNumber(), 0.0);
  double bestLL = -std::numeric_limits<double>::infinity();
  for (auto node: tree.getNodes()) {
    auto id = getSubtreeID(node);
    if (idToLL.find(id) != idToLL.end()) {
      // we have a likelihood value
      auto value = idToLL[id];
      nodeIdToLL[node->node_index] = value;
      bestLL = std::max<double>(value, bestLL);
    }
  }
  for (auto node: tree.getNodes()) {
    bool hasValue = nodeIdToLL[node->node_index] != 0.0;
    std::string label;
    if (hasValue) {
      double value = nodeIdToLL[node->node_index] - bestLL;
      if(!node->left && node->label) {
        label = std::string(node->label);
        free(node->label);
        node->label = nullptr;
        label += "_";
      }
      label += std::to_string(value);
      node->label = (char*)malloc(label.size() + 1);
      memcpy(node->label, label.c_str(), label.size());
      node->label[label.size()] = 0;
    }
  }
}



static void rootSearchAux(SpeciesTree &speciesTree, 
    SpeciesTreeLikelihoodEvaluator &evaluator,
    std::vector<unsigned int> &movesHistory, 
    std::vector<unsigned int> &bestMovesHistory, 
    double &bestLL, 
    unsigned int &visits,
    unsigned int maxDepth, 
    RootLikelihoods *rootLikelihoods,
    TreePerFamLLVec *treePerFamLLVec) 
{
  if (movesHistory.size() > maxDepth) {
    return;
  }
  std::vector<unsigned int> moves;
  moves.push_back(movesHistory.back() % 2);
  moves.push_back(2 + (movesHistory.back() % 2));
  for (auto direction: moves) {
    if (!SpeciesTreeOperator::canChangeRoot(speciesTree, direction)) {
      continue;
    }
    movesHistory.push_back(direction);
    evaluator.pushRollback();
    SpeciesTreeOperator::changeRoot(speciesTree, direction);
    evaluator.forceGeneRootOptimization();
    double ll = evaluator.computeLikelihood();
    if (treePerFamLLVec) {
      auto newick = speciesTree.getTree().getNewickString();
      treePerFamLLVec->push_back({newick, PerFamLL()});
      auto &perFamLL = treePerFamLLVec->back().second;
      evaluator.fillPerFamilyLikelihoods(perFamLL);
    }
    auto root = speciesTree.getRoot();
    if (rootLikelihoods) {
      rootLikelihoods->saveValue(root->left, ll);
      rootLikelihoods->saveValue(root->right, ll);
    }
    visits++;
    unsigned int additionalDepth = 0;
    if (ll > bestLL) {
      bestLL = ll;
      bestMovesHistory = movesHistory; 
      Logger::info << "Found better root " << ll << std::endl;
      additionalDepth = 3;
    }
    rootSearchAux(speciesTree, 
        evaluator,
        movesHistory, 
        bestMovesHistory, 
        bestLL, 
        visits,
        maxDepth + additionalDepth,
        rootLikelihoods,
        treePerFamLLVec);
    SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
    evaluator.popAndApplyRollback();
    movesHistory.pop_back();
  }
}

double SpeciesRootSearch::rootSearch(
    SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluator &evaluator,
    unsigned int maxDepth,
    RootLikelihoods *rootLikelihoods,
    TreePerFamLLVec *treePerFamLLVec)
{
  Logger::timed << "[Species search] Root search with depth=" << maxDepth << std::endl;
  std::vector<unsigned int> movesHistory;
  std::vector<unsigned int> bestMovesHistory;
  double bestLL = evaluator.computeLikelihood();
  if (treePerFamLLVec) {
    treePerFamLLVec->clear();
    auto newick = speciesTree.getTree().getNewickString();
    treePerFamLLVec->push_back({newick, PerFamLL()});
    auto &perFamLL = treePerFamLLVec->back().second;
    evaluator.fillPerFamilyLikelihoods(perFamLL);
  }
  auto root = speciesTree.getRoot();
  if (rootLikelihoods) {
    rootLikelihoods->saveValue(root->left, bestLL);
    rootLikelihoods->saveValue(root->right, bestLL);
  }
  unsigned int visits = 1;
  movesHistory.push_back(1);
  rootSearchAux(speciesTree,
      evaluator,
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits, 
      maxDepth,
      rootLikelihoods,
      treePerFamLLVec);
  movesHistory[0] = 0;
  rootSearchAux(speciesTree, 
      evaluator,
      movesHistory, 
      bestMovesHistory, 
      bestLL, 
      visits,
      maxDepth,
      rootLikelihoods,
      treePerFamLLVec);
  for (unsigned int i = 1; i < bestMovesHistory.size(); ++i) {
    SpeciesTreeOperator::changeRoot(speciesTree, bestMovesHistory[i]);
  }
  evaluator.forceGeneRootOptimization();
  {
    auto newick = speciesTree.getTree().getNewickString();
    PLLRootedTree tree(newick, false); 
    rootLikelihoods->fillTree(tree);
  }
  Logger::timed << "[Species search] After root search: LL=" << bestLL << std::endl;
  return bestLL;
}






