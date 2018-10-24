#ifndef _MOVES_H_
#define _MOVES_H_

#include <likelihoods/LibpllEvaluation.hpp>
#include <treeSearch/Rollbacks.hpp>
#include <memory>

class JointTree;

class Move {
public:
  virtual std::shared_ptr<Rollback> applyMove(JointTree &tree) = 0;
  static std::shared_ptr<Move> createNNIMove(int nodeIndex, bool left, bool blo);
  static std::shared_ptr<Move> createSPRMove(int pruneIndex, int regraftIndex, const vector<int> &path);
    
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
  SPRMove(int pruneIndex, int regraftIndex, const vector<int> &path);
  virtual std::shared_ptr<Rollback> applyMove(JointTree &tree);
  virtual std::ostream& print(std::ostream & os) const;
private:
  int pruneIndex_;
  int regraftIndex_;
  vector<int> path_;
};

#endif
