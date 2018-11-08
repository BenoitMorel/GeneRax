#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>
#include <chrono>
#include <ParallelContext.hpp>
#include "Arguments.hpp"
#include<limits>

int getTreeHashRec(pll_unode_t *node, int depth = 1) {
  int res = node->node_index * depth;
  if (node->next) {
    res += getTreeHashRec(node->next->back, depth + 1);
    res += getTreeHashRec(node->next->next->back, depth + 2);
  }
  return res;
}
    
int JointTree::getTreeHash()
{
  return getTreeHashRec(getTreeInfo()->root);
}

void printLibpllNode(pll_unode_s *node, Logger &os, bool isRoot)
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

void printLibpllTreeRooted(pll_unode_t *root, Logger &os){
  os << "(";
  printLibpllNode(root, os, true);
  os << ",";
  printLibpllNode(root->back, os, true);
  os << ");" << endl;
}


JointTree::JointTree(const string &newick_file,
    const string &alignment_file,
    const string &speciestree_file,
    const string &geneSpeciesMap_file,
    double dupRate,
    double lossRate):
  geneSpeciesMap_(geneSpeciesMap_file),
  dupRate_(dupRate),
  lossRate_(lossRate)
{
   info_.alignmentFilename = alignment_file;
  info_.model = "GTR";
  libpllEvaluation_ = LibpllEvaluation::buildFromFile(newick_file, info_);
  pllSpeciesTree_ = pll_rtree_parse_newick(speciestree_file.c_str());
  assert(pllSpeciesTree_);
  reconciliationEvaluation_ = make_shared<ReconciliationEvaluation>(pllSpeciesTree_,  geneSpeciesMap_);
  setRates(dupRate, lossRate);

}

void JointTree::printLibpllTree() const {
  printLibpllTreeRooted(libpllEvaluation_->getTreeInfo()->root, Logger::info);
}



void JointTree::optimizeParameters() {
  libpllEvaluation_->optimizeAllParameters();
  if (Arguments::costsEstimation) {
    optimizeDTRates();
  }
}

double JointTree::computeLibpllLoglk() {
  return libpllEvaluation_->computeLikelihood();
}

double JointTree::computeReconciliationLoglk () {
  return reconciliationEvaluation_->evaluate(libpllEvaluation_->getTreeInfo());
}

double JointTree::computeJointLoglk() {
  return computeLibpllLoglk() + computeReconciliationLoglk();
}

void JointTree::printLoglk(bool libpll, bool rec, bool joint, Logger &os) {
  if (joint)
    os << "joint: " << computeJointLoglk() << "  ";
  if (libpll)
    os << "libpll: " << computeLibpllLoglk() << "  ";
  if (rec)
    os << "reconciliation: " << computeReconciliationLoglk() << "  ";
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
  ofstream os(fileName);
  char *newick = pll_utree_export_newick(getTreeInfo()->root, 0);
  os << newick;
}

shared_ptr<pllmod_treeinfo_t> JointTree::getTreeInfo() {
  return libpllEvaluation_->getTreeInfo();
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
    double newLL = computeReconciliationLoglk();
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
  reconciliationEvaluation_->setRates(dup, loss);
}

