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

JointTree::JointTree(const std::string &newickString,
    const std::string &alignmentFilename,
    const std::string &speciestree_file,
    const std::string &geneSpeciesMapfile,
    const std::string &substitutionModel,
    RecModel reconciliationModel,
    RecOpt reconciliationOpt,
    bool rootedGeneTree,
    double supportThreshold,
    double recWeight,
    bool safeMode,
    bool optimizeDTLRates,
    const Parameters &ratesVector):
  _libpllEvaluation(newickString, false, alignmentFilename, substitutionModel),
  _speciesTree(speciestree_file, true),
  _optimizeDTLRates(optimizeDTLRates),
  _safeMode(safeMode),
  _enableReconciliation(true),
  _enableLibpll(true),
  _recOpt(reconciliationOpt),
  _recWeight(recWeight),
  _supportThreshold(supportThreshold)
{

  _geneSpeciesMap.fill(geneSpeciesMapfile, newickString);
  reconciliationEvaluation_ = std::make_unique<ReconciliationEvaluation>(_speciesTree,  
      _geneSpeciesMap, 
      reconciliationModel,
      rootedGeneTree);
  setRates(ratesVector);

}

JointTree::~JointTree()
{
}


void JointTree::printLibpllTree() const {
  printLibpllTreeRooted(_libpllEvaluation.getTreeInfo()->root, Logger::info);
}



void JointTree::optimizeParameters(bool felsenstein, bool reconciliation) {
  if (felsenstein && _enableLibpll) {
    _libpllEvaluation.optimizeAllParameters();
  }
  if (reconciliation && _enableReconciliation && _optimizeDTLRates) {
    std::cerr << "optimize parametres per family" << std::endl;
    if (reconciliationEvaluation_->implementsTransfers()) {  
      PerFamilyDTLOptimizer::optimizeDTLRates(*this, _recOpt);
    } else {
      PerFamilyDTLOptimizer::optimizeDLRates(*this, _recOpt);
    }
  }
}

double JointTree::computeLibpllLoglk(bool incremental) {
  if (!_enableLibpll) {
    return 1.0;
  }
  return _libpllEvaluation.computeLikelihood(incremental);
}

double JointTree::computeReconciliationLoglk () {
  if (!_enableReconciliation) {
    return 1.0;
  }
  return reconciliationEvaluation_->evaluate(getGeneTree()) * _recWeight;
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


void JointTree::applyMove(Move &move) {
  _rollbacks.push(std::move(move.applyMove(*this)));
}

void JointTree::optimizeMove(Move &move) {
  if (_enableLibpll) {
    move.optimizeMove(*this);
  }
}


void JointTree::rollbackLastMove() {
  assert(!_rollbacks.empty());
  _rollbacks.top()->applyRollback();
  _rollbacks.pop();
}

void JointTree::save(const std::string &fileName, bool append) {
  auto root = reconciliationEvaluation_->getRoot();
  if (!root) {
    root = reconciliationEvaluation_->inferMLRoot(getGeneTree());
  }
  assert(root);
  LibpllParsers::saveUtree(root, fileName, append);
}

pllmod_treeinfo_t * JointTree::getTreeInfo() {
  return _libpllEvaluation.getTreeInfo();
}


void JointTree::invalidateCLV(pll_unode_s *node)
{
  reconciliationEvaluation_->invalidateCLV(node->node_index);
  _libpllEvaluation.invalidateCLV(node->node_index);
}


void JointTree::setRates(const Parameters &ratesVector)
{
  _ratesVector = ratesVector;
  if (_enableReconciliation) {
    reconciliationEvaluation_->setRates(ratesVector);
  }
}


void JointTree::printInfo() 
{
  auto treeInfo = getTreeInfo();
  auto speciesLeaves = getSpeciesTree().getLeavesNumber();
  auto geneLeaves = treeInfo->tip_count;
  auto sites = treeInfo->partitions[0]->sites;
  Logger::info << "Species leaves: " << speciesLeaves << std::endl;
  Logger::info << "Gene leaves: " << geneLeaves << std::endl;
  Logger::info << "Sites: " << sites << std::endl;
  Logger::info << std::endl;
}

