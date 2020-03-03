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
class NeighborJoining {
public:
  /**
   * Infer a NJ tree, using the gene tree internode distances to compute
   * the distance matrix. The distance matrix is very similar to the one
   * built in NJst (another NJ took).
   */
  static std::unique_ptr<PLLRootedTree> geneTreeNJ(const Families &families);

  /**
   *  Infer a NJ tree, using the number of genes per species and per family
   *  to compute the distance matrix.
   *
   *  Warning: this method is not the best choice, and is very sensible to
   *  missing data.
   */
  static std::unique_ptr<PLLRootedTree> countProfileNJ(const Families &families);  
};
