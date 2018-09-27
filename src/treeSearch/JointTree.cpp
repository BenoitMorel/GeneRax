#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>
#include <chrono>


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
    double dupCost,
    double lossCost):
  dupCost_(dupCost),
  lossCost_(lossCost),
  transferCost_(0.01)
{

  info_.alignmentFilename = alignment_file;
  info_.model = "GTR";
  evaluation_ = LibpllEvaluation::buildFromFile(newick_file, info_);
  updateBPPTree();
  speciesTree_ = IO::readTreeFile(speciestree_file);
  vector<BPPTree> geneTrees(1, geneTree_);
  map_ = SpeciesGeneMapper::map(
      geneTrees.begin(), geneTrees.end(), *speciesTree_, trees)[0];
}

JointTree::JointTree(BPPTree geneTree,
    const LibpllAlignmentInfo *alignment,
    BPPTree speciesTree,
    const SpeciesGeneMap &map,
    double dupCost,
    double lossCost):
  evaluation_(LibpllEvaluation::buildFromPhylo(geneTree, *alignment)),
  speciesTree_(speciesTree),
  map_(map),
  info_(*alignment),
  dupCost_(dupCost),
  lossCost_(lossCost),
  transferCost_(0.0)
{
  updateBPPTree();
}


void JointTree::printLibpllTree() const {
  printLibpllTreeRooted(evaluation_->getTreeInfo()->root, cout);
}

void JointTree::printBPPTree() const {
  IO::write(*geneTree_, cout);
}

void JointTree::printSpeciesTree() const {
  IO::write(*speciesTree_, cout);
}

void JointTree::optimizeParameters() {
  evaluation_->optimizeAllParameters();
}

double JointTree::computeLibpllLoglk() {
  return evaluation_->computeLikelihood();
}

double JointTree::computeALELoglk () {
  auto genetree_copy = PhyloTreeToolBox::cloneTree(*geneTree_);
  PhyloTreeToolBox::removeArtificialGenes(*genetree_copy);
  double ale_loglk = ALEevaluation::evaluate(*genetree_copy, *speciesTree_, map_, 1, 1, dupCost_, transferCost_, lossCost_);
  return ale_loglk;
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
  ofstream os(fileName);
  IO::write(*geneTree_, os);
}

void JointTree::updateBPPTree() {
  geneTree_ = buildFromLibpll(evaluation_, evaluation_->getTreeInfo()->root);
}

BPPTree JointTree::getGeneTree() {
  return geneTree_;
}

shared_ptr<pllmod_treeinfo_t> JointTree::getTreeInfo() {
  return evaluation_->getTreeInfo();
}

JointTree& JointTree::getThreadInstance() {
  return *this;
}

ParallelJointTree::ParallelJointTree(BPPTree geneTree,
    const LibpllAlignmentInfo *alignment,
    BPPTree speciesTree,
    const SpeciesGeneMap &map,
    double dupCost,
    double lossCost,
    int threads)
{
  for (int i = 0; i < threads; ++i) {
    trees_.push_back(make_shared<JointTree>(geneTree,
          alignment, speciesTree, map, dupCost, lossCost));
  }
}

ParallelJointTree::ParallelJointTree(const string &newick_file,
    const string &alignment_file,
    const string &speciestree_file,
    double dupCost,
    double lossCost,
    int threads)
{
  for (int i = 0; i < threads; ++i) {
    trees_.push_back(make_shared<JointTree>(newick_file,
          alignment_file, speciestree_file,dupCost, lossCost));
  }
}


void ParallelJointTree::optimizeParameters() {
#pragma omp parallel for num_threads(getThreadsNumber())
  for (int i = 0; i < trees_.size(); ++i) {
    trees_[i]->optimizeParameters();
  }
}

int ParallelJointTree::getThreadsNumber() const {
  return trees_.size();
}

JointTree& ParallelJointTree::getThreadInstance() {
  int tid = omp_get_thread_num();
  if (tid >= trees_.size()) {
    cerr << "invalid index " << trees_.size() << " in getThreadInstance" << endl;
    exit(1);
  }
  return *trees_[tid];
}

void ParallelJointTree::applyMove(shared_ptr<Move> move) {
  for (auto &tree: trees_) {
    tree->applyMove(move);
  }
}

bool ParallelJointTree::checkConsistency() {
  cerr << "check consistency" << endl;
  vector<double> ll(getThreadsNumber());
#pragma omp parallel for num_threads(getThreadsNumber())
  for (int i = 0; i < trees_.size(); ++i) {
    ll[i] = trees_[i]->computeJointLoglk();
  }
  auto refll = ll[0];
  for (int i = 0; i < trees_.size(); ++i) {
    if (ll[i] != refll) {
      cerr << "Error, one tree at least has a different ll" << endl;
      exit(1);
    }
  }
}




