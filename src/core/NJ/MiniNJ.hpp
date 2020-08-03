#pragma once

#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>
#include <unordered_map>
#include <vector>
#include <string>
#include <util/types.hpp>

class PLLUnrootedTree;
class PLLRootedTree;
class GeneSpeciesMapping;

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


  static void computeDistanceMatrix(const Families &families,
      bool minMode, 
      bool reweight,
      bool useBL,
      bool useBootstrap,
      bool ustar,
      DistanceMatrix &distanceMatrix,
      std::vector<std::string> &speciesIdToSpeciesString,
      StringToUint &speciesStringToSpeciesId);

private:
  static std::unique_ptr<PLLRootedTree> geneTreeNJ(const Families &families, bool minAlgo, bool ustarAlgo = false, bool reweight = false);
};
