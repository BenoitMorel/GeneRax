#include "Rollbacks.hpp"
#include <treeSearch/JointTree.h>


void NNIRollback::applyRollback() {
  assert(PLL_SUCCESS == pllmod_tree_rollback(&rollback_));
  if (resetRoot_) {
    tree_.getReconciliationEvaluation()->setRoot(0);
  }
}

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
  auto prune = (pll_unode_s*)(rollback_.SPR.prune_edge);
  auto regraft = (pll_unode_s*)(rollback_.SPR.regraft_edge);
  tree_.invalidateCLV(prune->next->back);
  tree_.invalidateCLV(prune->next->next);
  tree_.invalidateCLV(prune->next);
  tree_.invalidateCLV(prune->next->next->back);
  tree_.invalidateCLV(regraft);
  tree_.invalidateCLV(regraft->back);
}

