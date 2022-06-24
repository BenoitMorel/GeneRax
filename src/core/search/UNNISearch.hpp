#pragma once
#include <trees/PLLUnrootedTree.hpp>


/**
 * Swap B and C
 */
struct UNNIMove {
  UNNIMove (corax_unode_t *e,
      bool l): edge(e), left(l) {}
  void apply();
  corax_unode_t *getA() const {return edge->next->back;}
  corax_unode_t *getB() const {return edge->next->next->back;}
  corax_unode_t *getC() const {return left? edge->back->next->back : 
    edge->back->next->next->back;}
  corax_unode_t *getD() const {return !left? edge->back->next->back : 
    edge->back->next->next->back;}

  corax_unode_t *edge;
  bool left;
};

/**
 *  Interface for the object evaluating
 *  the score to be optimized for an unrooted
 *  tree search.
 */
class USearchEvaluator {
public:
  virtual ~USearchEvaluator() {}
  virtual double eval(PLLUnrootedTree &tree) = 0;
  virtual double evalNNI(PLLUnrootedTree &tree, UNNIMove &move) = 0;
};


/**
 *
 */
class UNNISearch {
public:
  UNNISearch(PLLUnrootedTree &tree,
      USearchEvaluator &evaluator):
    _tree(tree),
    _evaluator(evaluator) {}
  double runRound();
  double search();
private:
  PLLUnrootedTree &_tree;
  USearchEvaluator &_evaluator;
};
