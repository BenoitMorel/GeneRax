#pragma once

#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>
#include <vector>
#include <string>
#include <util/types.hpp>

class NeighborJoining {
public:
  /**
   *  Apply Neighbor Joining algorithm
   *  @param DistanceMatrix a symetric distance matrix
   *  @param speciesIdToSpeciesString mapping from indices in 
   *    the distance matrix to the species labels
   *  @param speciesStringToSpeciesId mapping from the species
   *    labels to the indices in the distance matrix
   *  @param constrainTree optional constrain tree: if set, the NJ
   *    is only used to compute the same tree but with branch lengths
   *  @return The NJ tree (rooted, although the root is irrelevant)
   *
   *  Some parameters are passed by value on purpose
   *  
   */
  static std::unique_ptr<PLLRootedTree> applyNJ(
      DistanceMatrix distanceMatrix,
      std::vector<std::string> speciesIdToSpeciesString,
      StringToUint speciesStringToSpeciesId,
      PLLRootedTree *constrainTree = nullptr);
};


