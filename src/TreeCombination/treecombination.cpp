

#include <string>
#include <fstream>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/ParallelOfstream.hpp>
#include <parallelization/ParallelContext.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <maths/bitvector.hpp>

#include <unordered_map>

static void optimizeParameters(LibpllEvaluation &evaluation, double radius) 
{
  evaluation.optimizeAllParameters(radius);
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


double evalQuick(const std::string &treePath,
    bool isFile,
    const std::string &alignmentFile,
    const std::string &model,
    std::unordered_set<size_t> &cache)
{
  PLLUnrootedTree tree(treePath, isFile);
  auto hash = tree.getUnrootedTreeHash();
  if (cache.find(hash) != cache.end()) {
    return -std::numeric_limits<double>::infinity();
  }
  cache.insert(hash);
  LibpllEvaluation evaluation(treePath, isFile, alignmentFile, model);
  double ll = 0.0;
  ll = evaluation.computeLikelihood();
  evaluation.optimizeBranches(10.0, 1.0);
  ll = evaluation.computeLikelihood();
  Logger::timed << "Approx ll=\t" << ll << std::endl;
  return ll;
}
double eval(const std::string &treePath,
    bool isFile,
    const std::string &alignmentFile,
    const std::string &model,
    std::unordered_set<size_t> &cache,
    double &bestLL,
    std::string &bestTree,
    std::string &bestModel)
{
  PLLUnrootedTree tree(treePath, isFile);
  auto hash = tree.getUnrootedTreeHash();
  if (cache.find(hash) != cache.end()) {
    return -std::numeric_limits<double>::infinity();
  }
  cache.insert(hash);
  LibpllEvaluation evaluation(treePath, isFile, alignmentFile, model);
  
  double ll = 0.0;
  ll = evaluation.computeLikelihood();
  evaluation.optimizeBranches(10.0, 1.0);
  ll = evaluation.computeLikelihood();
  optimizeParameters(evaluation, 0.1);
  ll = evaluation.computeLikelihood();
  Logger::timed << "Eval ll= \t" << ll << std::endl;
  if (ll > bestLL) {
    Logger::timed << "Found better tree! ll=" << ll << std::endl;
    bestLL = ll;
    bestTree = evaluation.getGeneTree().getNewickString();
    bestModel = evaluation.getModelStr();
  }
  return ll;
}

double evalThorough(const std::string &treePath,
    bool isFile,
    const std::string &alignmentFile,
    const std::string &model,
    bool topology = false,
    std::string out = "")
{
  LibpllEvaluation evaluation(treePath, isFile, alignmentFile, model);
  optimizeParameters(evaluation, 0.1);
  evaluation.optimizeBranches(10.0, 1.0);
  if (topology) {
    optimizeTopology(evaluation, 5, 1, 20, 0.1);
  }
  if (out.size()) {
    evaluation.getGeneTree().save(out);
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
  
  std::string newick = "(";
  newick += subtree1;
  newick += ",";
  newick += subtree2;
  newick += ");";
  return newick;
}

void getBranchesToRecombine(PLLUnrootedTree &tree1,
    SplitToBranch &splitToBranch1,
    SplitToBranch &splitToBranch2,
    BranchToSplit &branchToSplit1,
    BranchPairs &branches)
{
  auto nodeNumber = tree1.getDirectedNodesNumber();
  auto nodes1 = tree1.getPostOrderNodes();
  std::vector<bool> isConflicting(nodeNumber, false);
  std::vector<bool> isSeriousConflicting(nodeNumber, false);
  std::vector<bool> hasConflictingChild(nodeNumber, false);
  for (auto node: nodes1) {
    if (!node->next || !node->back->next) {
      continue;
    }
    auto sp1 = branchToSplit1[node->node_index];
    if (splitToBranch2.find(sp1) == splitToBranch2.end()) {
      isConflicting[node->node_index] = true;
      if (node->length > 0.00001) {
        isSeriousConflicting[node->node_index] = true;
      }
    }
    auto leftIndex = node->next->back->node_index;
    auto rightIndex = node->next->next->back->node_index;
    if (isSeriousConflicting[leftIndex] || isSeriousConflicting[rightIndex]
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
    Logger::info << "lightsearch trees alignment model outputtree" 
      << std::endl;
    return 1;
  }
  int i = 1;
  std::string inputTreesPath(argv[i++]);
  std::string alignmentFile(argv[i++]);
  std::string model(argv[i++]);
  std::string outputTreePath(argv[i++]);
 
  std::ifstream is(inputTreesPath);
  std::vector<std::string> inputTreeStrings;
  std::string line;
  while(std::getline(is, line)) {
    inputTreeStrings.push_back(line);
  }
 
  std::cout.precision(11);

  assert(inputTreeStrings.size());
  PLLUnrootedTree tree1(inputTreeStrings[0], false);

  LabelToIndex labelToIndex;
  for (auto leaf: tree1.getLeaves()) {
    std::string label = leaf->label;
    unsigned int index = labelToIndex.size();
    labelToIndex.insert({label, index});
  }


  std::string bestTreeStr = inputTreeStrings[0];
  std::string bestModelStr = model;
  std::unordered_set<size_t> cache;
  std::unordered_set<size_t> cacheApprox;
  double bestLL = -99999999999999;
  Logger::info << "Evaluating all input trees..." << std::endl;
  for (const auto &treeStr: inputTreeStrings) {
    eval(treeStr, false, alignmentFile, bestModelStr, cacheApprox, bestLL, bestTreeStr, bestModelStr); 
  }
  auto initialLL = bestLL;
  for (unsigned int i = 1; i < inputTreeStrings.size(); ++i) {
    Logger::info << "Starting new round " << i << std::endl;
    const auto &newTreePath = inputTreeStrings[i];
    SplitToBranch splitToBranch1;
    SplitToBranch splitToBranch2;
    BranchToSplit branchToSplit1;
    BranchToSplit branchToSplit2;
    PLLUnrootedTree bestTree(bestTreeStr, false); 
    PLLUnrootedTree newTree(newTreePath, false);
    computeSplits(bestTree, splitToBranch1, branchToSplit1, labelToIndex);
    computeSplits(newTree, splitToBranch2, branchToSplit2, labelToIndex);
    BranchPairs branches;
    getBranchesToRecombine(bestTree, splitToBranch1,
        splitToBranch2, branchToSplit1, branches);
    
    
    
    for (auto pair: branches) {
      for (auto reverse: {true, false}) {
        auto recombinedTree = getRecombinedString(pair, reverse);
        evalQuick(recombinedTree, false, alignmentFile, bestModelStr, cacheApprox);
        eval(recombinedTree, false, alignmentFile, model, cache,
            bestLL, bestTreeStr, bestModelStr);
      }
    }
  }

  
  auto finalLL = evalThorough(bestTreeStr, false, alignmentFile, model, false, outputTreePath);
  Logger::info << std::endl;
  Logger::timed << "Initial LL " << initialLL << std::endl;
  Logger::timed << "Final LL " << finalLL << std::endl;
  
  ParallelContext::finalize();
  return 0;
}

int main(int argc, char** argv) {
  return lightSearch(argc, argv, nullptr);
}
