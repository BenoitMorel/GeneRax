#pragma once

#include <likelihoods/LibpllEvaluation.hpp>

#include <memory>
#include <vector>



class JointTree;


class SavedBranch {
public:
  SavedBranch(corax_unode_s *branch): 
    branch_(branch), 
    length_(branch->length) 
  {}

  void restore();
  corax_unode_t *getNode() {return branch_;}
private:
  corax_unode_t *branch_;
  double length_;
};

class SPRRollback {
public:
  virtual ~SPRRollback() {}
  SPRRollback(JointTree &tree, 
      corax_tree_rollback_t &rollback,
      const std::vector<SavedBranch> &branches,
      corax_unode_t *root):
    tree_(tree),
    rollback_(rollback),
    branches_(branches),
    root_(root)
  {}
  
  void saveBranch(const SavedBranch &branch)
  {
    branches_.push_back(branch);
  }

  virtual void applyRollback();

private:
  JointTree &tree_;
  corax_tree_rollback_t rollback_;
  std::vector<SavedBranch> branches_;
  corax_unode_t *root_;
};

