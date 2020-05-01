#pragma once

#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>

/*
 * Naive NJ implementation. 
 *
 * Both distance matrix and NJ tree computations could be 
 * implemented more efficiently if needed.
 */
class MiniNJ {
public:
  /**
   * Infer a NJ tree, using the gene tree internode distances to compute
   * the distance matrix. The distance matrix is very similar to the one
   * built in NJst (another NJ took).
   */
  static std::unique_ptr<PLLRootedTree> runMiniNJ(const Families &families);
  
  
  /**
   *  Run the original NJst algorithm
   */
  static std::unique_ptr<PLLRootedTree> runNJst(const Families &families);
  static std::unique_ptr<PLLRootedTree> runWeighted(const Families &families);
  
private:
  static std::unique_ptr<PLLRootedTree> geneTreeNJ(const Families &families, bool minAlgo, bool weightedAlgo);
};
