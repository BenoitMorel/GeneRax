#pragma once

#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>
#include <unordered_map>
#include <vector>
#include <string>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <algorithm>
#include <memory>
#include <util/types.hpp>


static const double invalidDouble = std::numeric_limits<double>::infinity();

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
  static std::unique_ptr<PLLRootedTree> runUstar(const Families &families);
  static std::unique_ptr<PLLRootedTree> runWMinNJ(const Families &families);
  
  static std::unique_ptr<PLLRootedTree> applyNJ(DistanceMatrix &distanceMatrix,
    std::vector<std::string> &speciesIdToSpeciesString,
    StringToUint &speciesStringToSpeciesId);
private:
  static std::unique_ptr<PLLRootedTree> geneTreeNJ(const Families &families, bool minAlgo, bool ustarAlgo = false, bool reweight = false);
};
