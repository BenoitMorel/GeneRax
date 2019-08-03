#include "JointTree.hpp"
#include <search/Moves.hpp>
#include <parallelization/ParallelContext.hpp>
#include <optimizers/PerFamilyDTLOptimizer.hpp>
#include <IO/LibpllParsers.hpp>
#include <chrono>
#include <limits>
#include <functional>


size_t leafHash(pll_unode_t *leaf) {
  assert(leaf);
  std::hash<std::string> hash_fn;
  return hash_fn(std::string(leaf->label));
}

size_t getTreeHashRec(pll_unode_t *node, size_t i) {
  assert(node);
  if (i == 0) 
    i = 1;
  if (!node->next) {
    return leafHash(node);
  }
  auto hash1 = getTreeHashRec(node->next->back, i + 1);
  auto hash2 = getTreeHashRec(node->next->next->back, i + 1);
  //Logger::info << "(" << hash1 << "," << hash2 << ") ";
  std::hash<size_t> hash_fn;
  auto m = std::min(hash1, hash2);
  auto M = std::max(hash1, hash2);
  return hash_fn(m * i + M);

}

pll_unode_t *findMinimumHashLeafRec(pll_unode_t * root, size_t &hashValue)
{
  assert(root);
  if (!root->next) {
    hashValue = leafHash(root);
    return root;
  }
  auto n1 = root->next->back;
  auto n2 = root->next->next->back;
  size_t hash1, hash2;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    hashValue = hash1;
    return min1;
  } else {
    hashValue = hash2;
    return min2;
  }
}

pll_unode_t *findMinimumHashLeaf(pll_unode_t * root) 
{
  assert(root);
  auto n1 = root;
  auto n2 = root->back;
  size_t hash1 = 0;
  size_t hash2 = 0;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    return min1;
  } else {
    return min2;
  }
}

void JointTree::printAllNodes(std::ostream &os)
{
  auto treeinfo = getTreeInfo();
  for (unsigned int i = 0; i < treeinfo->subnode_count; ++i) {
    auto node = treeinfo->subnodes[i];
    os << "node:" << node->node_index << " back:" << node->back->node_index;
    if (node->next) {
      os << " left:" << node->next->node_index << " right:" << node->next->next->node_index  << std::endl;
    } else {
      os << " label:" << node->label << std::endl;
    }
  }
}
    
size_t JointTree::getUnrootedTreeHash()
{
  auto minHashLeaf = findMinimumHashLeaf(getTreeInfo()->root);
  auto res = getTreeHashRec(minHashLeaf, 0) + getTreeHashRec(minHashLeaf->back, 0);
  return res % 100000;
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
  os << ");" << std::endl;
}


JointTree::JointTree(const std::string &newick_string,
    const std::string &alignment_file,
    const std::string &speciestree_file,
    const std::string &geneSpeciesMap_file,
    const std::string &substitutionModel,
    RecModel reconciliationModel,
    RecOpt reconciliationOpt,
    bool rootedGeneTree,
    double recWeight,
    bool safeMode,
    bool optimizeDTLRates,
    const DTLRatesVector &ratesVector):
  optimizeDTLRates_(optimizeDTLRates),
  safeMode_(safeMode),
  enableReconciliation_(true),
  enableLibpll_(true),
  recOpt_(reconciliationOpt),
  _recWeight(recWeight)
{
   info_.alignmentFilename = alignment_file;
  info_.model = substitutionModel;
  libpllEvaluation_ = LibpllEvaluation::buildFromString(newick_string, info_.alignmentFilename, info_.model);
  pllSpeciesTree_ = pll_rtree_parse_newick(speciestree_file.c_str());
  if (!pllSpeciesTree_) {
    Logger::perrank << "Cannot read tree " << speciestree_file << " !" << std::endl;
  }
  assert(pllSpeciesTree_);
  geneSpeciesMap_.fill(geneSpeciesMap_file, newick_string);
  reconciliationEvaluation_ = std::make_shared<ReconciliationEvaluation>(pllSpeciesTree_,  
      geneSpeciesMap_, 
      reconciliationModel,
      rootedGeneTree);
  setRates(ratesVector);

}

JointTree::~JointTree()
{
  pll_rtree_destroy(pllSpeciesTree_, 0);
  pllSpeciesTree_ = 0;
}


void JointTree::printLibpllTree() const {
  printLibpllTreeRooted(libpllEvaluation_->getTreeInfo()->root, Logger::info);
}



void JointTree::optimizeParameters(bool felsenstein, bool reconciliation) {
  if (felsenstein && enableLibpll_) {
    libpllEvaluation_->optimizeAllParameters();
  }
  if (reconciliation && enableReconciliation_ && optimizeDTLRates_) {
    if (reconciliationEvaluation_->implementsTransfers()) {  
      PerFamilyDTLOptimizer::optimizeDTLRates(*this, recOpt_);
    } else {
      PerFamilyDTLOptimizer::optimizeDLRates(*this, recOpt_);
    }
  }
}

double JointTree::computeLibpllLoglk(bool incremental) {
  if (!enableLibpll_) {
    return 1.0;
  }
  return libpllEvaluation_->computeLikelihood(incremental);
}

double JointTree::computeReconciliationLoglk () {
  if (!enableReconciliation_) {
    return 1.0;
  }
  return reconciliationEvaluation_->evaluate(libpllEvaluation_->getTreeInfo()) * _recWeight;
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
  os << std::endl;
}


pll_unode_t *JointTree::getNode(unsigned int index) {
  return getTreeInfo()->subnodes[index];
}


void JointTree::applyMove(std::shared_ptr<Move> move) {
  auto rollback = move->applyMove(*this);
  rollbacks_.push(rollback);
}

void JointTree::optimizeMove(std::shared_ptr<Move> move) {
  if (enableLibpll_) {
    move->optimizeMove(*this);
  }
}


void JointTree::rollbackLastMove() {
  assert(!rollbacks_.empty());
  rollbacks_.top()->applyRollback();
  rollbacks_.pop();
}

void JointTree::save(const std::string &fileName, bool append) {
  LibpllParsers::saveUtree(reconciliationEvaluation_->getRoot(), fileName, append);
}

std::shared_ptr<pllmod_treeinfo_t> JointTree::getTreeInfo() {
  return libpllEvaluation_->getTreeInfo();
}


void JointTree::invalidateCLV(pll_unode_s *node)
{
  reconciliationEvaluation_->invalidateCLV(node->node_index);
  libpllEvaluation_->invalidateCLV(node->node_index);
}


void JointTree::setRates(const DTLRatesVector &ratesVector)
{
  _ratesVector = ratesVector;
  if (enableReconciliation_) {
    reconciliationEvaluation_->setRates(ratesVector);
  }
}


void JointTree::printInfo() 
{
  auto treeInfo = getTreeInfo();
  auto speciesLeaves = getSpeciesTree()->tip_count;
  auto geneLeaves = treeInfo->tip_count;;
  auto sites = treeInfo->partitions[0]->sites;
  Logger::info << "Species leaves: " << speciesLeaves << std::endl;
  Logger::info << "Gene leaves: " << geneLeaves << std::endl;
  Logger::info << "Sites: " << sites << std::endl;
  Logger::info << std::endl;
}

