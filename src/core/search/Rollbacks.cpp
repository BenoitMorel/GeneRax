#include <search/Rollbacks.hpp>
#include <trees/JointTree.hpp>


void SavedBranch::restore() {
  pllmod_utree_set_length(branch_, length_);
}

void SPRRollback::applyRollback() {
  assert(PLL_SUCCESS == pllmod_tree_rollback(&rollback_));
  for (auto &b: branches_) {
    b.restore();
    tree_.invalidateCLV(b.getNode());
    tree_.invalidateCLV(b.getNode()->back);
  }
  auto prune = static_cast<pll_unode_s *>(rollback_.SPR.prune_edge);
  auto regraft = static_cast<pll_unode_s*>(rollback_.SPR.regraft_edge);
  tree_.invalidateCLV(prune->next->back);
  tree_.invalidateCLV(prune->next->next);
  tree_.invalidateCLV(prune->next);
  tree_.invalidateCLV(prune->next->next->back);
  tree_.invalidateCLV(regraft);
  tree_.invalidateCLV(regraft->back);
  tree_.setRoot(root_);
}

