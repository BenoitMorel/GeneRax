#include "Rollbacks.hpp"
#include <treeSearch/JointTree.h>


void NNIRollback::applyRollback() {
  assert(PLL_SUCCESS == pllmod_tree_rollback(&rollback_));
  tree_.updateBPPTree();
}

void SPRRollback::applyRollback() {
  assert(PLL_SUCCESS == pllmod_tree_rollback(&rollback_));
  tree_.updateBPPTree();
}

