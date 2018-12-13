#ifndef _ROLLBACKS_H_
#define _ROLLBACKS_H_

#include <likelihoods/LibpllEvaluation.hpp>
#include <memory>
#include <vector>

using namespace std;

class JointTree;

class Rollback {
public:
  virtual void applyRollback() = 0;
};

class NNIRollback: public Rollback {
public:
  NNIRollback(JointTree &tree, pll_tree_rollback_t &rollback, bool resetRoot):
    tree_(tree),
    rollback_(rollback),
    resetRoot_(resetRoot)
  {}

  virtual void applyRollback();

private:
  JointTree &tree_;
  pll_tree_rollback_t rollback_;
  bool resetRoot_;
};

class SavedBranch {
public:
  SavedBranch(pll_unode_s *branch): 
    branch_(branch), 
    length_(branch->length) 
  {}

  void restore();
  pll_unode_t *getNode() {return branch_;}
private:
  pll_unode_t *branch_;
  double length_;
};

class SPRRollback: public Rollback {
public:
  SPRRollback(JointTree &tree, 
      pll_tree_rollback_t &rollback,
      const vector<SavedBranch> &branches):
    tree_(tree),
    rollback_(rollback),
    branches_(branches)
  {}
  
  virtual void applyRollback();

private:
  JointTree &tree_;
  pll_tree_rollback_t rollback_;
  vector<SavedBranch> branches_;
};


#endif
