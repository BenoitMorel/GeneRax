/**
 * This is a simpler reimplementation of the raxml-ng search algorithm with autodetect, fast and slow SPR roudns
 */


#include <string>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/ParallelOfstream.hpp>
#include <parallelization/ParallelContext.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <maths/bitvector.hpp>

#include <unordered_map>

static void optimizeParameters(LibpllEvaluation &evaluation, double radius) 
{
  double initialLL = evaluation.computeLikelihood(false);
  //Logger::timed << "[" << initialLL << "]" << " Optimize parametres (" << radius << ")" << std::endl;
  
  evaluation.optimizeAllParameters(radius);
}

static void optimizeBranches(LibpllEvaluation &evaluation, double radius) 
{
  evaluation.optimizeBranches(radius);
}


static bool optimizeTopology(LibpllEvaluation &evaluation, 
    unsigned int radius, 
    unsigned int thorough, 
    unsigned int toKeep, 
    double cutoff) 
{
  
  double initialLL = evaluation.computeLikelihood(false);
  Logger::timed << "["  << initialLL << "] " 
    << (thorough ? "SLOW" : "FAST") << " SPR Round (radius=" << radius << ")" << std::endl;
  double ll = evaluation.raxmlSPRRounds(1, radius, thorough, toKeep, cutoff);
  Logger::info << "ll=" << ll << std::endl;
  return ll > initialLL + 0.1;
}
double eval(const std::string &treePath,
    bool isFile,
    const std::string &alignmentFile,
    const std::string &model,
    bool topology = false)
{
  LibpllEvaluation evaluation(treePath, isFile, alignmentFile, model);
  //optimizeParameters(evaluation, 0.1);
  optimizeBranches(evaluation, 0.1);
  Logger::timed << "ll=" << evaluation.computeLikelihood() << std::endl;
  
  if (topology) {
    optimizeParameters(evaluation, 0.1);
    optimizeTopology(evaluation, 25, 1, 20, 0.1);
  }
  return evaluation.computeLikelihood();
}

double evalModel(const std::string &treePath,
    bool isFile,
    const std::string &alignmentFile,
    std::string &model)
{
  LibpllEvaluation evaluation(treePath, isFile, alignmentFile, model);
  optimizeParameters(evaluation, 0.1);
  model = evaluation.getModelStr();
  Logger::timed << "ll=" << evaluation.computeLikelihood() << std::endl;
  return evaluation.computeLikelihood();
}



using Split = genesis::utils::Bitvector;
using LabelToIndex = std::unordered_map<std::string, unsigned int>;
using SplitToBranch = std::unordered_map<Split, corax_unode_t *>;
using BranchToSplit = std::unordered_map<unsigned int, Split>;
using BranchPair = std::pair<corax_unode_t *, corax_unode_t *>;
using BranchPairs = std::vector<BranchPair>;


Split computeSplitsRec(corax_unode_t *branch,
    SplitToBranch &splitToBranch,
    LabelToIndex &labelToIndex)
{
  auto N = labelToIndex.size();
  if (!branch->next) {
    // leaf case
    Split split(N, false);
    split.set(labelToIndex[branch->label]);
    return split;
  }
  // internal branch case
  auto left = branch->next->back;
  auto right = branch->next->next->back;
  auto sp1 = computeSplitsRec(left, splitToBranch, labelToIndex);
  auto sp2 = computeSplitsRec(right, splitToBranch, labelToIndex);
  auto split = sp1 | sp2;
  if (split[0]) {
    splitToBranch.insert({split, branch});
  } else {
    splitToBranch.insert({~split, branch->back});
  }
  return split;
}

void computeSplits(PLLUnrootedTree &tree,
    SplitToBranch &splitToBranch,
    BranchToSplit &branchToSplit,
    LabelToIndex &labelToIndex)
{
  // find the reference branch
  corax_unode_t *refLeaf = nullptr;
  for (auto leaf: tree.getLeaves()) {
    if (labelToIndex[leaf->label] == 0) {
      refLeaf = leaf;
      break;
    }
  }
  assert(refLeaf);
  computeSplitsRec(refLeaf->back, splitToBranch, labelToIndex);
  for (auto it: splitToBranch) {
    branchToSplit.insert({it.second->node_index, it.first});
    branchToSplit.insert({it.second->back->node_index, it.first});
  }
}

std::string getRecombinedString(const BranchPair &pair,
    bool reverse)
{
  auto n1 = pair.first;
  auto n2 = pair.second;
  if (reverse) {
    std::swap(n1, n2);
  }
  
  std::string subtree1 = PLLUnrootedTree::getSubtreeString(n1);
  std::string subtree2 = PLLUnrootedTree::getSubtreeString(n2->back);
 
  /*
  if (subtree1.find("ENSSTUP00000053904") != std::string::npos) {
    subtree1 = PLLUnrootedTree::getSubtreeString(n1->back);
  }
  if (subtree2.find("ENSSTUP00000053904") == std::string::npos) {
    subtree2 = PLLUnrootedTree::getSubtreeString(n2->back);
  }
  */
  std::string newick = "(";
  newick += subtree1;
  newick += ",";
  newick += subtree2;
  newick += ");";
  return newick;
}

void getBranchesToRecombine(PLLUnrootedTree &tree1,
    PLLUnrootedTree &tree2,
    SplitToBranch &splitToBranch1,
    SplitToBranch &splitToBranch2,
    BranchToSplit &branchToSplit1,
    BranchToSplit &branchToSplit2,
    BranchPairs &branches)
{
  auto nodeNumber = tree1.getDirectedNodesNumber();
  auto nodes1 = tree1.getPostOrderNodes();
  std::vector<bool> isConflicting(nodeNumber, false);
  std::vector<bool> hasConflictingChild(nodeNumber, false);
  for (auto node: nodes1) {
    if (!node->next || !node->back->next) {
      continue;
    }
    auto sp1 = branchToSplit1[node->node_index];
    if (splitToBranch2.find(sp1) == splitToBranch2.end()) {
      isConflicting[node->node_index] = true;
    }
    auto leftIndex = node->next->back->node_index;
    auto rightIndex = node->next->next->back->node_index;
    if (isConflicting[leftIndex] || isConflicting[rightIndex]
        || hasConflictingChild[leftIndex]
        || hasConflictingChild[rightIndex]) {
      hasConflictingChild[node->node_index] = true;
    }
  }
  for (auto branch: tree1.getBranches()) {    
    auto index = branch->node_index;
    auto indexBack = branch->back->node_index;
    assert(isConflicting[index] == isConflicting[indexBack]);
    if (isConflicting[index]) {
      continue;
    }
    if (hasConflictingChild[index] && hasConflictingChild[indexBack]) {
      auto sp = branchToSplit1[index];
      auto branch1 = splitToBranch1[sp];
      auto branch2 = splitToBranch2[sp];
      branches.push_back(BranchPair(branch1, branch2));
    }
  }
}




int lightSearch(int argc, char** argv, void* comm)
{
  ParallelContext::init(comm);
  Logger::init();
  if (argc != 5) {
    Logger::info << "Syntax error" << std::endl;
    Logger::info << "lightsearch tree1 tree2 alignment model" 
      << std::endl;
    return 1;
  }
  int i = 1;
  std::string treePath1(argv[i++]);
  std::string treePath2(argv[i++]);
  std::string alignmentFile(argv[i++]);
  std::string model(argv[i++]);
  
  

  PLLUnrootedTree tree1(treePath1);
  PLLUnrootedTree tree2(treePath2);

  LabelToIndex labelToIndex;
  for (auto leaf: tree1.getLeaves()) {
    std::string label = leaf->label;
    unsigned int index = labelToIndex.size();
    labelToIndex.insert({label, index});
  }
  
  SplitToBranch splitToBranch1;
  SplitToBranch splitToBranch2;
  BranchToSplit branchToSplit1;
  BranchToSplit branchToSplit2;

  computeSplits(tree1, splitToBranch1, branchToSplit1, labelToIndex);
  computeSplits(tree2, splitToBranch2, branchToSplit2, labelToIndex);
 
  Logger::info << "Number of leaves: " << tree1.getLeavesNumber() << std::endl;
  unsigned int RF = 0;
  for (auto it: splitToBranch1) {
    RF += splitToBranch2.find(it.first) == splitToBranch2.end();
  }
  Logger::info << "RF=" << RF << std::endl;
 

  
  std::string fixedModel = model;
  evalModel(treePath1, true, alignmentFile, fixedModel);  
  eval(treePath2, true, alignmentFile, fixedModel);  
  
  BranchPairs branches;
  getBranchesToRecombine(tree1, tree2, splitToBranch1,
      splitToBranch2, branchToSplit1, branchToSplit2, branches);

  Logger::info << "Branches: " << branches.size() << std::endl;

  double bestLL = -99999999999999;
  std::string bestTree;
  for (auto pair: branches) {
    for (auto reverse: {true}) {// {true, false}) {
      auto newTreeStr1 = getRecombinedString(pair, reverse);
      auto ll1 = eval(newTreeStr1, false, alignmentFile, fixedModel);
      if (ll1 > bestLL) {
        bestLL = ll1;
        bestTree = newTreeStr1;
      }
    }
  }
  //Logger::info << bestTree << std::endl;
  Logger::info << " best ll = " << bestLL << std::endl;
  
  
  auto finalLL = eval(bestTree, false, alignmentFile, model, true);
  Logger::info << "Final LL " << finalLL << std::endl;
  ParallelContext::finalize();
  return 0;
}

int main(int argc, char** argv) {
  return lightSearch(argc, argv, nullptr);
}
