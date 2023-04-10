#include "SpeciesSearchCommon.hpp"

#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>

static std::string getSubtreeID(corax_rnode_t *subtree)
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

void RootLikelihoods::saveValue(corax_rnode_t *subtree, double ll) 
{
  auto id = getSubtreeID(subtree); 
  idToLL[id] = ll;
}

void RootLikelihoods::fillTree(PLLRootedTree &tree)
{
  std::vector<double> nodeIdToLL(tree.getNodeNumber(), 0.0);
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

void SpeciesSearchState::betterTreeCallback(double ll)
{
  bool masterRankOnly = true;
  speciesTree.saveToFile(pathToBestSpeciesTree, masterRankOnly);
  bestLL = ll;
}

bool SpeciesSearchCommon::testSPR(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int prune,
    unsigned int regraft
    )
{
  evaluation.pushRollback();
  // Apply the move
  auto rollback = SpeciesTreeOperator::applySPRMove(speciesTree, prune, regraft);
  bool runExactTest = true;
  double approxLL = 0.0;
  if (evaluation.providesFastLikelihoodImpl()) {
    // first test with approximative likelihood
    approxLL = evaluation.computeLikelihoodFast();
    if (searchState.averageApproxError.isSignificant()) {
      //  Decide whether we can already
      // discard the move
      auto epsilon = 2.0 * searchState.averageApproxError.getAverage();
      runExactTest &= (approxLL + epsilon > searchState.bestLL );
    } 
  }
  if (runExactTest) {
    // we test the move with exact likelihood
    double testedTreeLL = evaluation.computeLikelihood();
    if (evaluation.providesFastLikelihoodImpl()) {
      searchState.averageApproxError.addValue(testedTreeLL - approxLL);
    }
    if (testedTreeLL > searchState.bestLL + 0.00000001) {
      searchState.betterTreeCallback(testedTreeLL);
      // Better tree found! Do not rollback, and return
      return true;
    }
  }
  // the tree is not better, rollback the move
  SpeciesTreeOperator::reverseSPRMove(speciesTree, prune, rollback);
  evaluation.popAndApplyRollback();
  return false;
}

bool SpeciesSearchCommon::veryLocalSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState,
    unsigned int spid)
{

  const unsigned int radius = 2;
  std::vector<unsigned int> prunes;
  prunes.push_back(spid);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(speciesTree, 
        prune, 
        radius, 
        regrafts);
    for (auto regraft: regrafts) {
      if (testSPR(speciesTree, evaluation, searchState, prune, regraft)) {
        Logger::timed << "\tfound better* (LL=" 
            << searchState.bestLL << ", hash=" << 
            speciesTree.getHash() << ")" << std::endl;
        veryLocalSearch(speciesTree, evaluation, searchState,
            prune);
        return true;
      }
    }
  }
  return false;
}



