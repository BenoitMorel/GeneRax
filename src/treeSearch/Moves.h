#ifndef _MOVES_H_
#define _MOVES_H_

#include <likelihoods/LibpllEvaluation.h>
#include <memory>

class JointTree;

class Rollback {
public:
  Rollback(JointTree &tree, pll_tree_rollback_t &rollback):
    tree_(tree),
    rollback_(rollback) 
  {}

  void applyRollback();

private:
  JointTree &tree_;
  pll_tree_rollback_t rollback_;

};


class Move {
public:
  virtual std::shared_ptr<Rollback> applyMove(JointTree &tree) = 0;
  static std::shared_ptr<Move> createNNIMove(int nodeIndex, bool left, bool blo);
  static std::shared_ptr<Move> createSPRMove(int pruneIndex, int regraftIndex);
    
  friend std::ostream & operator <<( std::ostream &os, const Move &move ) {
    return move.print(os);
  }
  
  virtual std::ostream& print(std::ostream & os) const = 0;
};

class NNIMove: public Move {
public:
    NNIMove(int nodeIndex, bool left, bool blo);
    virtual std::shared_ptr<Rollback> applyMove(JointTree &tree);
    virtual std::ostream& print(std::ostream & os) const;

private:
    //Optimizes the 5 branches around the NNI move
    void applyBLO(JointTree &tree, pll_unode_s *edge);
    int nodeIndex_;
    bool left_;
    bool blo_;
};

class SPRMove: public Move {
public:
  SPRMove(int pruneIndex, int regraftIndex);
  virtual std::shared_ptr<Rollback> applyMove(JointTree &tree);
  virtual std::ostream& print(std::ostream & os) const;
private:
  void applyBLO(JointTree &tree, 
    pll_unode_t *prune,
    pll_unode_t *regraft);
  int pruneIndex_;
  int regraftIndex_;
};

#endif
