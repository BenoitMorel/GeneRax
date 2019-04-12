#include <treeSearch/Moves.hpp>
#include <trees/JointTree.hpp>
#include <IO/Logger.hpp>

// constants taken from RAXML
#define DEF_LH_EPSILON            0.1
#define OPT_LH_EPSILON            0.1
#define RAXML_PARAM_EPSILON       0.001  //0.01
#define RAXML_BFGS_FACTOR         1e7
#define RAXML_BRLEN_SMOOTHINGS    5
#define RAXML_BRLEN_DEFAULT       0.1
#define RAXML_BRLEN_MIN           1.0e-6
#define RAXML_BRLEN_MAX           100.
#define RAXML_BRLEN_TOLERANCE     1.0e-9
#define RAXML_FREERATE_MIN        0.001
#define RAXML_FREERATE_MAX        100.
#define RAXML_BRLEN_SCALER_MIN    0.01
#define RAXML_BRLEN_SCALER_MAX 100.
  


shared_ptr<Move> Move::createSPRMove(int pruneIndex, int regraftIndex, const vector<int> &path) {
  return make_shared<SPRMove>(pruneIndex, regraftIndex, path);
}


void optimizeBranchesSlow(JointTree &tree,
    const vector<pll_unode_t *> &nodesToOptimize)
{
    auto root = tree.getTreeInfo()->root;
    // could be incremental and thus faster
    unsigned int params_indices[4] = {0,0,0,0};
    auto treeinfo = tree.getTreeInfo();
    
    tree.computeLibpllLoglk(); // update CLVs
    for (int j = 0; j < 2; ++j) 
    for (unsigned int i = 0; i < nodesToOptimize.size(); ++i) {
        pllmod_treeinfo_set_root(treeinfo.get(), nodesToOptimize[i]);
        double oldLoglk = tree.computeLibpllLoglk(true);
        double newLoglk = pllmod_opt_optimize_branch_lengths_local(
            treeinfo->partitions[0],
            treeinfo->root,
            params_indices,
            RAXML_BRLEN_MIN,
            RAXML_BRLEN_MAX,
            RAXML_BRLEN_TOLERANCE,
            RAXML_BRLEN_SMOOTHINGS,
            0,
            true);
       assert(oldLoglk <= newLoglk);
    }
    pllmod_treeinfo_set_root(treeinfo.get(), root);

}

bool equals(pll_unode_t *node1, pll_unode_t *node2) {
  return node1 == node2 || (node1->next && (node1->next == node2 || node1->next->next == node2));
}

SPRMove::SPRMove(int pruneIndex, int regraftIndex, const vector<int> &path):
  pruneIndex_(pruneIndex),
  regraftIndex_(regraftIndex),
  path_(path)
{
}

void printNode(pll_unode_s *node)
{
  Logger::info << node;
  if (node->next) {
    Logger::info << " " << node->next << " " << node->next->next;
  }
  Logger::info << endl;
}



shared_ptr<Rollback> SPRMove::applyMove(JointTree &tree)
{
  auto root = tree.getRoot();
  auto prune = tree.getNode(pruneIndex_);
  auto regraft = tree.getNode(regraftIndex_);
  tree.invalidateCLV(prune->next->back);
  tree.invalidateCLV(prune->next->next);
  tree.invalidateCLV(prune->next);
  tree.invalidateCLV(prune->next->next->back);
  tree.invalidateCLV(regraft);
  tree.invalidateCLV(regraft->back);
  for (int branchIndex: path_) {
    tree.invalidateCLV(tree.getNode(branchIndex));
    tree.invalidateCLV(tree.getNode(branchIndex)->back);
  }
  pll_tree_rollback_t pll_rollback;
  vector<SavedBranch> savedBranches;;
    
  branchesToOptimize_.push_back(prune);
  branchesToOptimize_.push_back(regraft->back);
  branchesToOptimize_.push_back(regraft);
  for (int branchIndex: path_) {
    auto node = tree.getNode(branchIndex);
    branchesToOptimize_.push_back(node);
    if (path_.size() == 1) {
      branchesToOptimize_.push_back(node->next);
      branchesToOptimize_.push_back(node->next->next);
      node = node->back;
      branchesToOptimize_.push_back(node->next);
      branchesToOptimize_.push_back(node->next->next);
    }
  }
  for (auto branch: branchesToOptimize_) {
    savedBranches.push_back(branch);
  }
  assert(PLL_SUCCESS == pllmod_utree_spr(prune, regraft, &pll_rollback));
  rollback_ = make_shared<SPRRollback>(tree, pll_rollback, savedBranches, root);
  return rollback_;
}
  
void SPRMove::optimizeMove(JointTree &tree)
{
  optimizeBranchesSlow(tree, branchesToOptimize_);
  branchesToOptimize_.clear();
}

ostream& SPRMove::print(ostream & os) const {
  os << "SPR(";
  os << "prune:" <<pruneIndex_ << ", regraft:" << regraftIndex_;
  os << ",path_size:" << path_.size();
  os << ")";
  return os;
}


