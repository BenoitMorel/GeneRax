#include "Rollbacks.hpp"
#include <treeSearch/JointTree.h>


void NNIRollback::applyRollback() {
  assert(PLL_SUCCESS == pllmod_tree_rollback(&rollback_));
  tree_.updateBPPTree();
}

void SavedBranch::restore() {
  pllmod_utree_set_length(branch_, length_);
}


void SPRRollback::applyRollback() {
  assert(PLL_SUCCESS == pllmod_tree_rollback(&rollback_));
  for (auto &b: branches_)
    b.restore();
  tree_.updateBPPTree();
}

