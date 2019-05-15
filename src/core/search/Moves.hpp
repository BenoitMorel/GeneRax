#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <search/Rollbacks.hpp>

#include <memory>



class JointTree;

class Move {
public:
  virtual ~Move() {}
  
  static std::shared_ptr<Move> createSPRMove(int pruneIndex, int regraftIndex, const std::vector<int> &path);
  
  virtual std::shared_ptr<Rollback> applyMove(JointTree &tree) = 0;
  
  virtual void optimizeMove(JointTree &tree) = 0;
    
  friend std::ostream & operator <<( std::ostream &os, const Move &move ) {
    return move.print(os);
  }
  
  virtual std::ostream& print(std::ostream & os) const = 0;
};

class SPRMove: public Move {
public:
  SPRMove(int pruneIndex, int regraftIndex, const std::vector<int> &path);
  virtual ~SPRMove() {}
  virtual std::shared_ptr<Rollback> applyMove(JointTree &tree);
  virtual void optimizeMove(JointTree &tree);
  virtual std::ostream& print(std::ostream & os) const;
private:
  int pruneIndex_;
  int regraftIndex_;
  std::vector<int> path_;
  std::vector<pll_unode_t *> branchesToOptimize_;
  std::shared_ptr<SPRRollback> rollback_;
};
