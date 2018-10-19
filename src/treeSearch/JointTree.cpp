#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>
#include <chrono>
#include <ParallelContext.hpp>
#include "Arguments.hpp"
#include<limits>

void printLibpllNode(pll_unode_s *node, ostream &os, bool isRoot)
{
  if (node->next) {
    os << "(";
    printLibpllNode(node->next->back, os, false);
    os << ",";
    printLibpllNode(node->next->next->back, os, false);
    os << ")";
  } else {
    os << node->label;
  }
  os << ":" << (isRoot ? node->length / 2.0 : node->length);
}

void printLibpllTreeRooted(pll_unode_t *root, ostream &os){
  os << "(";
  printLibpllNode(root, os, true);
  os << ",";
  printLibpllNode(root->back, os, true);
  os << ");" << endl;
}

BPPNode createNode(BPPTree tree, BPPNode father = 0, const std::string&name = "", double branchLength = 0.0) {
  BPPNode node(new bpp::PhyloNode(name));
  if (!father) {
    tree->createNode(node);
  } else {
    BPPBranch branch(new bpp::PhyloBranch);
    if (tree->hasFather(father)) {
      branch->setLength(branchLength);
    } else {
      // special case: this branch corresponds to half the branch
      // of the root in the unrooted tree
      branch->setLength(branchLength / 2.0);
    }
    tree->createNode(father, node, branch);
  }
  return node;
}

void addFromLibpll(BPPTree tree, BPPNode bppFather, pll_unode_s *libpllNode)
{
  BPPNode newNode = createNode(tree, bppFather, libpllNode->label ? libpllNode->label : "", libpllNode->length);
  if (libpllNode->next) {
    addFromLibpll(tree, newNode, libpllNode->next->back);
    addFromLibpll(tree, newNode, libpllNode->next->next->back);
  }
}

std::shared_ptr<bpp::PhyloTree> buildFromLibpll(std::shared_ptr<LibpllEvaluation> evaluation, pll_unode_s *libpllRoot)
{
  std::shared_ptr<bpp::PhyloTree> tree(new bpp::PhyloTree(true));
  BPPNode root = createNode(tree, 0, "root");
  tree->setRoot(root);
  addFromLibpll(tree, root, libpllRoot);
  addFromLibpll(tree, root, libpllRoot->back);
  return tree;
}

JointTree::JointTree(const string &newick_file,
    const string &alignment_file,
    const string &speciestree_file,
    double dupRate,
    double lossRate):
  dupRate_(dupRate),
  lossRate_(lossRate),
  transferRate_(0.0),
  aleWeight_(Arguments::aleWeight),
  needAleRecomputation_(true)
{
   info_.alignmentFilename = alignment_file;
  info_.model = "GTR";
  evaluation_ = LibpllEvaluation::buildFromFile(newick_file, info_);
  updateBPPTree();
  speciesTree_ = IO::readTreeFile(speciestree_file);
  vector<BPPTree> geneTrees(1, geneTree_);
  pllSpeciesTree_ = pll_rtree_parse_newick(speciestree_file.c_str());
  assert(pllSpeciesTree_);
  map_ = SpeciesGeneMapper::map(
      geneTrees.begin(), geneTrees.end(), *speciesTree_, trees)[0];
  aleEvaluation_ = make_shared<ALEEvaluation>(*speciesTree_, pllSpeciesTree_,  map_);
  setRates(dupRate, lossRate);

}

void JointTree::printLibpllTree() const {
  printLibpllTreeRooted(evaluation_->getTreeInfo()->root, Logger::info);
}

void JointTree::printBPPTree() const {
  IO::write(*geneTree_, Logger::info);
}

void JointTree::printSpeciesTree() const {
  IO::write(*speciesTree_, Logger::info);
}

void JointTree::optimizeParameters() {
  evaluation_->optimizeAllParameters();
}

double JointTree::computeLibpllLoglk() {
  return evaluation_->computeLikelihood();
}

double JointTree::computeALELoglk () {
  if (needAleRecomputation_) {
    aleLL_ = aleEvaluation_->evaluate(evaluation_->getTreeInfo());
    //aleLL_ = aleEvaluation_->evaluate(*geneTree_);
    needAleRecomputation_ = false;
  }
  return aleWeight_ * aleLL_;
}

double JointTree::computeJointLoglk() {
  return computeLibpllLoglk() + computeALELoglk();
}

void JointTree::printLoglk(bool libpll, bool ale, bool joint, ostream &os) {
  if (joint)
    os << "joint: " << computeJointLoglk() << "  ";
  if (libpll)
    os << "libpll: " << computeLibpllLoglk() << "  ";
  if (ale)
    os << "ale: " << computeALELoglk() << "  ";
  os << endl;
}


// todobenoit make it faster
pll_unode_t *JointTree::getNode(int index) {
  return getTreeInfo()->subnodes[index];
}


void JointTree::applyMove(shared_ptr<Move> move) {

  auto rollback = move->applyMove(*this);
  rollbacks_.push(rollback);
}


void JointTree::rollbackLastMove() {
  assert(!rollbacks_.empty());
  rollbacks_.top()->applyRollback();
  rollbacks_.pop();
}

void JointTree::save(const string &fileName) {
  updateBPPTree();
  ofstream os(fileName);
  IO::write(*geneTree_, os);
}

void JointTree::updateBPPTree() {
  geneTree_ = buildFromLibpll(evaluation_, evaluation_->getTreeInfo()->root);
  needAleRecomputation_ = true;
}

BPPTree JointTree::getGeneTree() {
  return geneTree_;
}

shared_ptr<pllmod_treeinfo_t> JointTree::getTreeInfo() {
  return evaluation_->getTreeInfo();
}

void JointTree::optimizeDTRates() {
  double bestLL = numeric_limits<double>::lowest();
  double bestDup = 0.0;
  double bestLoss = 0.0;
  double min = 0.001;
  double max = 2.0;
  int steps = 15;
  int begin = ParallelContext::getBegin(steps * steps);
  int end = ParallelContext::getEnd(steps * steps);
  
  for (int s = begin; s < end; ++s) {
    int i = s / 15;
    int j = s % 15;
    double dup = min + (max - min) * double(i) / double(steps);
    double loss = min + (max - min) * double(j) / double(steps);
    setRates(dup, loss);
    double newLL = computeALELoglk();
    if (newLL > bestLL) { 
      bestDup = dup;
      bestLoss = loss;
      bestLL = newLL;
    }
  }
  int bestRank = 0;
  ParallelContext::getBestLL(bestLL, bestRank);
  ParallelContext::broadcoastDouble(bestRank, bestDup);
  ParallelContext::broadcoastDouble(bestRank, bestLoss);
  setRates(bestDup, bestLoss);
  Logger::info << " best rates: " << bestDup << " " << bestLoss << endl;
}

void JointTree::setRates(double dup, double loss) { 
  dupRate_ = dup; 
  lossRate_ = loss;
  aleEvaluation_->setRates(dup, loss);
  needAleRecomputation_ = true;
}

