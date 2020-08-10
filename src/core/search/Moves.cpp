#include <search/Moves.hpp>
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
  


std::unique_ptr<Move> Move::createSPRMove(unsigned int pruneIndex, 
    unsigned int regraftIndex, 
    const std::vector<unsigned int> &path) {
  return std::make_unique<SPRMove>(pruneIndex, regraftIndex, path);
}


static void optimizeBranchesSlow(JointTree &tree,
    const std::vector<pll_unode_t *> &nodesToOptimize)
{
    auto root = tree.getTreeInfo()->root;
    // could be incremental and thus faster
    auto treeinfo = tree.getTreeInfo();
    auto ratecats = tree.getModel().num_ratecats();
    std::vector<unsigned int> params_indices(ratecats, 0);
    tree.computeLibpllLoglk(); // update CLVs
    for (unsigned int j = 0; j < 2; ++j) {
      for (unsigned int i = 0; i < nodesToOptimize.size(); ++i) {
          pllmod_treeinfo_set_root(treeinfo, nodesToOptimize[i]);
          double oldLoglk = tree.computeLibpllLoglk(true);
          double newLoglk = pllmod_opt_optimize_branch_lengths_local(
              treeinfo->partitions[0],
              treeinfo->root,
              &params_indices[0],
              RAXML_BRLEN_MIN,
              RAXML_BRLEN_MAX,
              RAXML_BRLEN_TOLERANCE,
              RAXML_BRLEN_SMOOTHINGS,
              0,
              true);
         assert(oldLoglk <= newLoglk);
      }
    }
    pllmod_treeinfo_set_root(treeinfo, root);
}


SPRMove::SPRMove(unsigned int pruneIndex, unsigned int regraftIndex, const std::vector<unsigned int> &path):
  pruneIndex_(pruneIndex),
  regraftIndex_(regraftIndex),
  path_(path)
{
}

std::unique_ptr<Rollback> SPRMove::applyMove(JointTree &tree)
{
  auto root = tree.getRoot();
  auto prune = tree.getNode(pruneIndex_);
  auto regraft = tree.getNode(regraftIndex_);
  assert (prune && prune->next);
  assert (regraft && regraft);
  tree.invalidateCLV(prune->next->back);
  tree.invalidateCLV(prune->next->next);
  tree.invalidateCLV(prune->next);
  tree.invalidateCLV(prune->next->next->back);
  tree.invalidateCLV(regraft);
  tree.invalidateCLV(regraft->back);
  for (auto branchIndex: path_) {
    tree.invalidateCLV(tree.getNode(branchIndex));
    tree.invalidateCLV(tree.getNode(branchIndex)->back);
  }
  pll_tree_rollback_t pll_rollback;
  std::vector<SavedBranch> savedBranches;
    
  branchesToOptimize_.push_back(prune);
  branchesToOptimize_.push_back(regraft->back);
  branchesToOptimize_.push_back(regraft);
  for (auto branchIndex: path_) {
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
  return std::make_unique<SPRRollback>(tree, pll_rollback, savedBranches, root);
}
  
void SPRMove::optimizeMove(JointTree &tree)
{
  optimizeBranchesSlow(tree, branchesToOptimize_);
  branchesToOptimize_.clear();
}

std::ostream& SPRMove::print(std::ostream & os) const {
  os << "SPR(";
  os << "prune:" <<pruneIndex_ << ", regraft:" << regraftIndex_;
  os << ",path_size:" << path_.size();
  os << ")";
  return os;
}


