#ifndef _ROLLBACKS_H_
#define _ROLLBACKS_H_

#include <likelihoods/LibpllEvaluation.h>
#include <memory>

class JointTree;

class Rollback {
public:
  virtual void applyRollback() = 0;
};

class NNIRollback: public Rollback {
public:
  NNIRollback(JointTree &tree, pll_tree_rollback_t &rollback):
    tree_(tree),
    rollback_(rollback) 
  {}

  virtual void applyRollback();

private:
  JointTree &tree_;
  pll_tree_rollback_t rollback_;
};


class SPRRollback: public Rollback {
public:
  SPRRollback(JointTree &tree, pll_tree_rollback_t &rollback):
    tree_(tree),
    rollback_(rollback) 
  {}

  virtual void applyRollback();

private:
  JointTree &tree_;
  pll_tree_rollback_t rollback_;
};


#endif
