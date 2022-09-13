

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
    std::unordered_map<size_t, double> &cache)
{
  PLLUnrootedTree tree(treePath, isFile);
  auto hash = tree.getUnrootedTreeHash();
  if (cache.find(hash) != cache.end()) {
    return cache[hash];
  }
  LibpllEvaluation evaluation(treePath, isFile, alignmentFile, model);
  double ll = 0.0;
  ll = evaluation.computeLikelihood();
  evaluation.optimizeBranches(10.0, 1.0);
  ll = evaluation.computeLikelihood();
  cache.insert({hash, ll});
  Logger::timed << "Approx ll=\t" << ll << std::endl;
  return ll;
}
double eval(const std::string &treePath,
    bool isFile,
    const std::string &alignmentFile,
    const std::string &model,
    std::unordered_map<size_t, double> &cache,
    double &bestLL,
    std::string &bestTree,
    std::string &bestModel)
{
  PLLUnrootedTree tree(treePath, isFile);
  auto hash = tree.getUnrootedTreeHash();
  if (cache.find(hash) != cache.end()) {
    return cache[hash];
  }
  LibpllEvaluation evaluation(treePath, isFile, alignmentFile, bestModel);
  Model m(model);
  evaluation.setParametersToOptimize(m.params_to_optimize());
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
    Logger::info << " best model " << bestModel << std::endl;
  }
  cache.insert({hash, ll});
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

int getNodesNumber(corax_unode_t *node)
{
  if (!node->next) {
    return 1;
  }
  return getNodesNumber(node->next->back) + getNodesNumber(node->next->next->back);
}


std::string getRecombinedString(const BranchPair &pair,
    bool reverse)
{
  auto n1 = pair.first;
  auto n2 = pair.second;
  if (reverse) {
    n1 = n1->back;
    n2 = n2->back;
  }
  //Logger::info << "n1: " << getNodesNumber(n1)  << "\t n2: " << getNodesNumber(n2->back) << std::endl;
  

  std::string subtree1 = PLLUnrootedTree::getSubtreeString(n1);
  std::string subtree2 = PLLUnrootedTree::getSubtreeString(n2->back);
  if (getNodesNumber(n1) < getNodesNumber(n1->back)) {
    return "";
  }
  std::string newick = "(";
  newick += subtree1;
  newick += ",";
  newick += subtree2;
  newick += ");";
  return newick;
}


std::string buildMultibinedSubtree(corax_unode_t *node2,
    const std::vector<corax_unode_t *> &branch2ToBranch1,
    bool &justPolytomies,
    bool firstCall = true) {
  if (!node2->next) {
    return std::string(node2->label);
  }
  if (!firstCall && branch2ToBranch1[node2->node_index]) {
    auto node1 = branch2ToBranch1[node2->node_index];
    return PLLUnrootedTree::getSubtreeString(node1);
  } else {
    if (node2->length > 0.0000011) {
      justPolytomies = false;
    }
    auto left = PLLUnrootedTree::getLeft(node2);
    auto right = PLLUnrootedTree::getRight(node2);
    auto leftString = buildMultibinedSubtree(left, branch2ToBranch1, justPolytomies, false);
    auto rightString = buildMultibinedSubtree(right, branch2ToBranch1, justPolytomies, false);
    std::string res = "(";
    res += leftString;
    res += ",";
    res += rightString;
    res += ")";
    return res;
  }
}


/*
 *  return directions such that directions[node->node_index] is
 *  true if node points to a node with label refLabel
 */
std::vector<bool> getBranchDirections(PLLUnrootedTree &tree,
    const std::string &refLabel)
{
  auto nodes = tree.getPostOrderNodes();
  std::vector<bool> directions(nodes.size(), false);
  for (auto node: nodes) {
    if (!node->next) {
      directions[node->node_index] = (refLabel == node->label);
    } else {
      auto leftIndex = PLLUnrootedTree::getLeft(node)->node_index;
      auto rightIndex = PLLUnrootedTree::getRight(node)->node_index;
      directions[node->node_index] = directions[leftIndex] ||
        directions[rightIndex];
    }
  }
  return directions;
}


std::vector<std::string> buildMultibinedTrees(PLLUnrootedTree &tree1,
    PLLUnrootedTree &tree2,
    BranchToSplit &branchToSplit1,
    SplitToBranch &splitToBranch1,
    SplitToBranch &splitToBranch2
  )

{
  std::vector<std::string> trees;

  auto refLabel = std::string(tree1.getAnyLeaf()->label);
  auto directions1 = getBranchDirections(tree1, refLabel);
  auto directions2 = getBranchDirections(tree2, refLabel);
  
  std::vector<corax_unode_t *> branch2ToBranch1(tree1.getDirectedNodesNumber(), nullptr);
  for (auto branch1: tree1.getPostOrderNodes()) {
    auto it = branchToSplit1.find(branch1->node_index);
    if (it == branchToSplit1.end()) {
      continue;
    }
    auto split = it->second;
    if (splitToBranch2.end() == splitToBranch2.find(split)) {
      continue;
    }
    auto branch2 = splitToBranch2.at(split);
    if (directions1[branch1->node_index] != directions2[branch2->node_index]) {
      branch2 = branch2->back;
    }
    branch2ToBranch1[branch2->node_index] = branch1;
  }

  for (auto branch2: tree2.getPostOrderNodes()) {
    if (!branch2->next || !branch2->back->next) {
      continue;
    }
    auto branch1 = branch2ToBranch1[branch2->node_index];
    if (!branch1) {
      continue;
    }
    // branch1 and branch2 agree and are not trivial branches
    bool justPolytomies = true;
    auto treeLeft = buildMultibinedSubtree(branch2, branch2ToBranch1, justPolytomies);
    if (justPolytomies) {
      continue;
    }
    auto treeRight = PLLUnrootedTree::getSubtreeString(branch1->back);
    std::string tree = "(";
    tree += treeLeft;
    tree += ",";
    tree += treeRight;
    tree += ");";
    //std::cout << "ADD TREE " << tree << std::endl;
    trees.push_back(tree);
  }
  Logger::info << "Added " << trees.size() << " trees" << std::endl;
  return trees;
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
  //ParallelContext::init(comm);
  Logger::init();
  if (argc != 5) {
    Logger::info << "Syntax error" << std::endl;
    Logger::info << "lightsearch trees alignment model prefix" 
      << std::endl;
    return 1;
  }
  int i = 1;
  std::string inputTreesPath(argv[i++]);
  std::string alignmentFile(argv[i++]);
  std::string model(argv[i++]);
  std::string prefix(argv[i++]);
  std::string outputLog = prefix + ".log";
  std::string outputStats = prefix + ".stats";
  std::string outputTreePath = prefix + ".newick";
  
  Logger::initFileOutput(outputLog);
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
  std::unordered_map<size_t, double> cache;
  std::unordered_map<size_t, double> cache2;
  std::unordered_map<size_t, double> cacheApprox;
  double bestLL = -99999999999999;
  Logger::info << "Evaluating all input trees..." << std::endl;
  unsigned int bestTreeNumber = 0;
  for (const auto &treeStr: inputTreeStrings) {
    auto previousBestLL = bestLL;
    const double EPS = 0.00001;
    auto ll = eval(treeStr, false, alignmentFile, model, cache, bestLL, bestTreeStr, bestModelStr);
    if (ll - previousBestLL > EPS) {
      bestTreeNumber = 1; // new best tree
    } else if (fabs(ll - previousBestLL) < EPS) {
      bestTreeNumber++; // this tree is as good as the best one
    }
  }
  Logger::info << "The best tree was found " << bestTreeNumber << " times." << std::endl;
 
  auto initialLL = bestLL;
  double finalLL = bestLL;
  if (bestTreeNumber < inputTreeStrings.size() / 2) {
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
#define MULTI
#define SIMPLE
#ifdef SIMPLE
      BranchPairs branches;
      getBranchesToRecombine(bestTree, splitToBranch1,
          splitToBranch2, branchToSplit1, branches);
      for (auto pair: branches) {
        for (auto reverse: {true, false}) {
          auto recombinedTree = getRecombinedString(pair, reverse);
          if (recombinedTree.size() == 0) {
            continue;
          }
          //evalQuick(recombinedTree, false, alignmentFile, bestModelStr, cacheApprox);
          eval(recombinedTree, false, alignmentFile, model, cache,
              bestLL, bestTreeStr, bestModelStr);
        }
      }
#endif
#ifdef MULTI
      bool improved = true;
      while (improved) {
        improved = false;
        PLLUnrootedTree currentBestTree(bestTreeStr, false);
        splitToBranch1.clear();
        branchToSplit1.clear();
        computeSplits(currentBestTree, splitToBranch1, branchToSplit1, labelToIndex);
        auto candidates = buildMultibinedTrees(currentBestTree,
          newTree,
          branchToSplit1,
          splitToBranch1,
          splitToBranch2
          );
        for (auto recombinedTree: candidates) {
          bool saveBestLL = bestLL;  
          //evalQuick(recombinedTree, false, alignmentFile, bestModelStr, cacheApprox);
          eval(recombinedTree, false, alignmentFile, model, cache,
                bestLL, bestTreeStr, bestModelStr);
          if (bestLL > saveBestLL) {
            improved = true;
            break;
          }
        }
      }
#endif
    }

  } else {  
    Logger::info << "Skipping the recombination step." << std::endl;
  }
  
  finalLL = evalThorough(bestTreeStr, false, alignmentFile, model, false, outputTreePath);
  Logger::info << std::endl;
  Logger::timed << "Initial LL " << initialLL << std::endl;
  Logger::timed << "Final LL " << finalLL << std::endl;
  std::ofstream statsOs(outputStats);
  statsOs << bestTreeNumber << " " << initialLL << " " << finalLL << std::endl; 
  //ParallelContext::finalize();
  return 0;
}

int main(int argc, char** argv) {
  return lightSearch(argc, argv, nullptr);
}
