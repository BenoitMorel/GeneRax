#include <treeSearch/JointTree.hpp>
#include <treeSearch/Moves.hpp>
#include <ParallelContext.hpp>

#include <chrono>
#include <limits>
#include <functional>


size_t leafHash(pll_unode_t *leaf) {
  hash<string> hash_fn;
  return hash_fn(string(leaf->label));
}

size_t getTreeHashRec(pll_unode_t *node, size_t i) {
  if (i == 0) 
    i = 1;
  if (!node->next) {
    return leafHash(node);
  }
  int hash1 = getTreeHashRec(node->next->back, i + 1);
  int hash2 = getTreeHashRec(node->next->next->back, i + 1);
  //Logger::info << "(" << hash1 << "," << hash2 << ") ";
  hash<int> hash_fn;
  int m = min(hash1, hash2);
  int M = max(hash1, hash2);
  return hash_fn(m * i + M);

}

pll_unode_t *findMinimumHashLeafRec(pll_unode_t * root, size_t &hash)
{
  if (!root->next) {
    hash = leafHash(root);
    return root;
  }
  auto n1 = root->next->back;
  auto n2 = root->next->next->back;
  size_t hash1, hash2;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    hash = hash1;
    return min1;
  } else {
    hash = hash2;
    return min2;
  }
}

pll_unode_t *findMinimumHashLeaf(pll_unode_t * root) 
{
  auto n1 = root;
  auto n2 = root->back;
  size_t hash1, hash2;
  auto min1 = findMinimumHashLeafRec(n1, hash1);
  auto min2 = findMinimumHashLeafRec(n2, hash2);
  if (hash1 < hash2) {
    return min1;
  } else {
    return min2;
  }
}

void JointTree::printAllNodes(ostream &os)
{
  auto treeinfo = getTreeInfo();
  for (unsigned int i = 0; i < treeinfo->subnode_count; ++i) {
    auto node = treeinfo->subnodes[i];
    os << "node:" << node->node_index << " back:" << node->back->node_index;
    if (node->next) {
      os << " left:" << node->next->node_index << " right:" << node->next->next->node_index  << endl;
    } else {
      os << " label:" << node->label << endl;
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
  os << ");" << endl;
}


JointTree::JointTree(const string &newick_string,
    const string &alignment_file,
    const string &speciestree_file,
    const string &geneSpeciesMap_file,
    const string &substitutionModel,
    const string &reconciliationModel,
    bool rootedGeneTree,
    bool safeMode,
    bool optimizeDTLRates,
    double dupRate,
    double lossRate,
    double transRate):
  geneSpeciesMap_(geneSpeciesMap_file),
  dupRate_(dupRate),
  lossRate_(lossRate),
  transRate_(transRate),
  optimizeDTLRates_(optimizeDTLRates),
  safeMode_(safeMode)
{
   info_.alignmentFilename = alignment_file;
  info_.model = substitutionModel;
  libpllEvaluation_ = LibpllEvaluation::buildFromString(newick_string, info_.alignmentFilename, info_.model);
  pllSpeciesTree_ = pll_rtree_parse_newick(speciestree_file.c_str());
  assert(pllSpeciesTree_);
  reconciliationEvaluation_ = make_shared<ReconciliationEvaluation>(pllSpeciesTree_,  
      geneSpeciesMap_, 
      reconciliationModel,
      rootedGeneTree);
  setRates(dupRate, lossRate, transRate);

}

void JointTree::printLibpllTree() const {
  printLibpllTreeRooted(libpllEvaluation_->getTreeInfo()->root, Logger::info);
}



void JointTree::optimizeParameters(bool felsenstein, bool reconciliation) {
  if (felsenstein) {
    libpllEvaluation_->optimizeAllParameters();
  }
  if (reconciliation) {
    if (reconciliationEvaluation_->implementsTransfers()) {  
      optimizeDTLRates();
    } else {
      optimizeDLRates();
    }
  }
}

double JointTree::computeLibpllLoglk(bool incremental) {
  return libpllEvaluation_->computeLikelihood(incremental);
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


pll_unode_t *JointTree::getNode(int index) {
  return getTreeInfo()->subnodes[index];
}


void JointTree::applyMove(shared_ptr<Move> move) {
  auto rollback = move->applyMove(*this);
  rollbacks_.push(rollback);
}

void JointTree::optimizeMove(shared_ptr<Move> move) {
  move->optimizeMove(*this);
}


void JointTree::rollbackLastMove() {
  assert(!rollbacks_.empty());
  rollbacks_.top()->applyRollback();
  rollbacks_.pop();
}

void JointTree::save(const string &fileName, bool append) {
  ofstream os(fileName, (append ? ofstream::app : ofstream::out));
  char *newick = pll_utree_export_newick(getTreeInfo()->root, 0);
  os << newick;
}

shared_ptr<pllmod_treeinfo_t> JointTree::getTreeInfo() {
  return libpllEvaluation_->getTreeInfo();
}

bool isValidLikelihood(double ll) {
  return isnormal(ll) && ll < 100000000000;
}

void JointTree::findBestRates(double minDup, double maxDup,
    double minLoss, double maxLoss, int steps,
    double &bestDup,
    double &bestLoss,
    double &bestLL) 
{
  bestLL = numeric_limits<double>::lowest();
  int totalSteps = pow(steps, 2);
  int begin = ParallelContext::getBegin(totalSteps);
  int end = ParallelContext::getEnd(totalSteps);
  for (int s = begin; s < end; ++s) {
    int i = s / steps;
    int j = s % steps;
    double dup = minDup + (maxDup - minDup) * double(i) / double(steps);
    double loss = minLoss + (maxLoss - minLoss) * double(j) / double(steps);
    setRates(dup, loss);
    double newLL = computeReconciliationLoglk();
    if (!isValidLikelihood(newLL)) {
      continue;
    }
    if (newLL > bestLL) { 
      bestDup = dup;
      bestLoss = loss;
      bestLL = newLL;
    }
  }
  int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  setRates(bestDup, bestLoss);
}

void JointTree::findBestRatesDTL(double minDup, double maxDup,
    double minLoss, double maxLoss, 
    double minTrans, double maxTrans, 
    int steps,
    double &bestDup,
    double &bestLoss,
    double &bestTrans,
    double &bestLL) 
{
  bestLL = numeric_limits<double>::lowest();
  int totalSteps = pow(steps, 3);
  int begin = ParallelContext::getBegin(totalSteps);
  int end = ParallelContext::getEnd(totalSteps);
  for (int s = begin; s < end; ++s) {
    int i = s / (steps * steps);
    int j = (s / steps) % steps;
    int k = s % steps;
    double dup = minDup + (maxDup - minDup) * double(i) / double(steps);
    double loss = minLoss + (maxLoss - minLoss) * double(j) / double(steps);
    double trans = minTrans + (maxTrans - minTrans) * double(k) / double(steps);
    setRates(dup, loss, trans);
    double newLL = computeReconciliationLoglk();
    if (!isValidLikelihood(newLL)) {
      continue;
    }
    if (newLL > bestLL) { 
      bestDup = dup;
      bestLoss = loss;
      bestTrans = trans;
      bestLL = newLL;
    }
  }
  int bestRank = 0;
  ParallelContext::getMax(bestLL, bestRank);
  ParallelContext::broadcastDouble(bestRank, bestDup);
  ParallelContext::broadcastDouble(bestRank, bestLoss);
  ParallelContext::broadcastDouble(bestRank, bestTrans);
  setRates(bestDup, bestLoss, bestTrans);
}


void JointTree::optimizeDLRates() {
  if (!optimizeDTLRates_) {
    Logger::info << "Skipping DL rates optimization (rates were defined by the user)" << endl;
    return;
  }
  Logger::timed << "Start optimizing DL rates" << endl;
  double bestLL = numeric_limits<double>::lowest();
  double newLL = 0;
  double bestDup = 0.0;
  double bestLoss = 0.0;
  double minDup = 0.0;
  double maxDup = 10.0;
  double minLoss = 0.0;
  double maxLoss = 10.0;
  int steps = 10;
  double epsilon = 0.001;
  do {
    bestLL = newLL;
    findBestRates(minDup, maxDup, minLoss, maxLoss, steps, bestDup, bestLoss, newLL);
    double offsetDup = (maxDup - minDup) / steps;
    double offsetLoss =(maxLoss - minLoss) / steps;
    minDup = max(0.0, bestDup - offsetDup);
    maxDup = bestDup + offsetDup;
    minLoss = max(0.0, bestLoss - offsetLoss);
    maxLoss = bestLoss + offsetLoss;
  } while (fabs(newLL - bestLL) > epsilon);
  Logger::info << " best rates: " << bestDup << " " << bestLoss <<  " " << newLL << endl;
  if  (!isValidLikelihood(newLL)) {
    Logger::error << "Invalid likelihood " << newLL << endl;
    ParallelContext::abort(10);
  }
}
    
void JointTree::optimizeDTLRates() {
  if (!optimizeDTLRates_) {
    Logger::info << "Skipping DL rates optimization (rates were defined by the user)" << endl;
    return;
  }
  Logger::timed << "Start optimizing DTL rates" << endl;
  double bestLL = numeric_limits<double>::lowest();
  double newLL = 0;
  double bestDup = 0.0;
  double bestLoss = 0.0;
  double bestTrans = 0.0;
  double minDup = 0.0;
  double maxDup = 1.0;
  double minLoss = 0.0;
  double maxLoss = 1.0;
  double minTrans = 0.0;
  double maxTrans = 1.0;
  int steps = 5;
  double epsilon = 0.01;
  do {
    bestLL = newLL;
    findBestRatesDTL(minDup, maxDup, minLoss, maxLoss, minTrans, maxTrans, steps, bestDup, bestLoss, bestTrans, newLL);
    double offsetDup = (maxDup - minDup) / steps;
    double offsetLoss = (maxLoss - minLoss) / steps;
    double offsetTrans = (maxTrans - minTrans) / steps;
    minDup = max(0.0, bestDup - offsetDup);
    maxDup = bestDup + offsetDup;
    minLoss = max(0.0, bestLoss - offsetLoss);
    maxLoss = bestLoss + offsetLoss;
    minTrans = max(0.0, bestTrans - offsetTrans);
    maxTrans = bestTrans + offsetTrans;
  } while (fabs(newLL - bestLL) > epsilon);
  if  (!isValidLikelihood(newLL)) {
    Logger::error << "Invalid likelihood " << newLL << endl;
    ParallelContext::abort(10);
  }
}
    
void JointTree::invalidateCLV(pll_unode_s *node)
{
  reconciliationEvaluation_->invalidateCLV(node->node_index);
  libpllEvaluation_->invalidateCLV(node->node_index);
}



void JointTree::setRates(double dup, double loss, double trans) { 
  dupRate_ = dup; 
  lossRate_ = loss;
  transRate_ = trans;
  reconciliationEvaluation_->setRates(dup, loss, trans);
}

void JointTree::printInfo() 
{
  auto treeInfo = getTreeInfo();
  int speciesLeaves = getSpeciesTree()->tip_count;
  int geneLeaves = treeInfo->tip_count;;
  int sites = treeInfo->partitions[0]->sites;
  Logger::info << "Species leaves: " << speciesLeaves << endl;
  Logger::info << "Gene leaves: " << geneLeaves << endl;
  Logger::info << "Sites: " << sites << endl;
  Logger::info << endl;
}

