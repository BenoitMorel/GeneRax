#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <treeSearch/Rollbacks.hpp>

#include <memory>

using namespace std;

class JointTree;

class Move {
public:
  static shared_ptr<Move> createSPRMove(int pruneIndex, int regraftIndex, const vector<int> &path);
  
  virtual shared_ptr<Rollback> applyMove(JointTree &tree) = 0;
  
  virtual void optimizeMove(JointTree &tree) = 0;
    
  friend ostream & operator <<( ostream &os, const Move &move ) {
    return move.print(os);
  }
  
  virtual ostream& print(ostream & os) const = 0;
};

class SPRMove: public Move {
public:
  SPRMove(int pruneIndex, int regraftIndex, const vector<int> &path);
  virtual shared_ptr<Rollback> applyMove(JointTree &tree);
  virtual void optimizeMove(JointTree &tree);
  virtual ostream& print(ostream & os) const;
private:
  int pruneIndex_;
  int regraftIndex_;
  vector<int> path_;
  vector<pll_unode_t *> branchesToOptimize_;
  shared_ptr<SPRRollback> rollback_;
};

