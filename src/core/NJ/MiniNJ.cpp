#include "MiniNJ.hpp"

#include "NeighborJoining.hpp"
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Families.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <algorithm>
#include <parallelization//ParallelContext.hpp>

void fillDistancesRec(pll_unode_t *currentNode, 
    bool useBL,
    bool useBootstrap,
    double currentDistance,
    std::vector<double> &distances)
{
  if (useBL) {
    currentDistance += std::max(0.0000001, currentNode->length);
  } else if (useBootstrap) {
    double bootstrapValue = (nullptr == currentNode->label) ? 0.0 : std::atof(currentNode->label);
    currentDistance += bootstrapValue;
  } else {
    // do not account for internal nodes which have too small
    // branch lengths
    if (currentNode->length == 0.0 || currentNode->length > 0.0000011) {
      currentDistance += 1.0;
    }
  }
  if (!currentNode->next) {
    // leaf
    distances[currentNode->node_index] = currentDistance;
    return;
  }
  fillDistancesRec(currentNode->next->back, useBL, useBootstrap, currentDistance, distances);
  fillDistancesRec(currentNode->next->next->back, useBL, useBootstrap, currentDistance, distances);
} 


static void geneDistancesFromGeneTree(PLLUnrootedTree &geneTree,
    GeneSpeciesMapping &mapping,
    StringToUint &speciesStringToSpeciesId,
    DistanceMatrix &distances,
    DistanceMatrix &distancesDenominator,
    bool minMode,
    bool reweight,
    bool useBL,
    bool useBootstrap,
    bool ustar)
{
  unsigned int speciesNumber = distances.size();
  if (speciesNumber < 3) {
    Logger::error << "Found less than 3 species, please check that the genes cover at least 3 species or that the gene-species mappings are correct" << std::endl;
    ParallelContext::abort(1);
  }
  std::vector<double>zeros(speciesNumber, 0.0);
  DistanceMatrix distancesToAdd(speciesNumber, zeros);
  DistanceMatrix distancesDenominatorToAdd(speciesNumber, zeros);

  auto leaves = geneTree.getLeaves();
  // build geneId -> speciesId
  std::vector<unsigned int> geneIdToSpeciesId(leaves.size());
  for (auto leafNode: leaves) {
    auto &species = mapping.getSpecies(leafNode->label);
    auto speciesId = speciesStringToSpeciesId[species];
    geneIdToSpeciesId[leafNode->node_index] = speciesId;
  }
  // build gene leaf distance matrix
  std::vector<double>zerosLeaf(leaves.size(), 0.0);
  std::vector<std::vector<double> > leafDistances(leaves.size(), zerosLeaf);
  for (auto leafNode: leaves) {
    fillDistancesRec(leafNode->back, useBL, useBootstrap, 0.0, leafDistances[leafNode->node_index]);
  }

  // fill species distance matrices
  for (auto gene1: leaves) {
    auto gid1 = gene1->node_index;
    auto spid1 = geneIdToSpeciesId[gid1];
    for (auto gene2: leaves) {
      auto gid2 = gene2->node_index;
      auto spid2 = geneIdToSpeciesId[gid2];
      if (!minMode) {
        distancesToAdd[spid1][spid2] += leafDistances[gid1][gid2];
        distancesDenominatorToAdd[spid1][spid2]++;
      } else {
        if (distancesDenominatorToAdd[spid1][spid2]) {
          distancesToAdd[spid1][spid2] =
            std::min(distancesToAdd[spid1][spid2], 
                   leafDistances[gid1][gid2]);
        } else {
          distancesToAdd[spid1][spid2] = leafDistances[gid1][gid2];
          distancesDenominatorToAdd[spid1][spid2] = 1;
        }
      }
    }
  }
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      if (reweight) {
        double factor = (double(leaves.size()));
        distancesToAdd[i][j] *= factor;
        distancesDenominatorToAdd[i][j] *= factor;
      }
      if (ustar) {
        if (distancesDenominatorToAdd[i][j] > 0.0) {
          distancesToAdd[i][j] /= distancesDenominatorToAdd[i][j];
          distancesDenominatorToAdd[i][j] = 1.0;
        }
      }
      distances[i][j] += distancesToAdd[i][j];
      distancesDenominator[i][j] += distancesDenominatorToAdd[i][j];
    }
  }
}

std::unique_ptr<PLLRootedTree> MiniNJ::runNJst(const Families &families)
{
  return geneTreeNJ(families, false);
}

std::unique_ptr<PLLRootedTree> MiniNJ::runWMinNJ(const Families &families)
{
  return geneTreeNJ(families, true, false, true);
}

std::unique_ptr<PLLRootedTree> MiniNJ::runUstar(const Families &families)
{
  return geneTreeNJ(families, false, true, false);
}

std::unique_ptr<PLLRootedTree> MiniNJ::runMiniNJ(const Families &families)
{
  return geneTreeNJ(families, true);
}


std::unique_ptr<PLLRootedTree> MiniNJ::geneTreeNJ(const Families &families, bool minMode, bool ustar, bool reweight)
{
  DistanceMatrix distanceMatrix;
  std::vector<std::string> speciesIdToSpeciesString;
  StringToUint speciesStringToSpeciesId;
  bool useBL = false;
  bool useBootstrap = false;
  computeDistanceMatrix(families,
      minMode,
      reweight,
      useBL,
      useBootstrap,
      ustar,
      distanceMatrix,
      speciesIdToSpeciesString,
      speciesStringToSpeciesId);
  auto res = NeighborJoining::applyNJ(distanceMatrix, 
      speciesIdToSpeciesString, 
      speciesStringToSpeciesId);
  /*
  return NeighborJoining::applyNJ(distanceMatrix, 
      speciesIdToSpeciesString, 
      speciesStringToSpeciesId,
      res.get());
    */
  return res;
}


void MiniNJ::computeDistanceMatrix(const Families &families,
  bool minMode, 
  bool reweight,
  bool useBL,
  bool useBootstrap,
  bool ustar,
  DistanceMatrix &distanceMatrix,
  std::vector<std::string> &speciesIdToSpeciesString,
  StringToUint &speciesStringToSpeciesId)
{
  for (auto &family: families) {
    GeneSpeciesMapping mappings;
    mappings.fill(family.mappingFile, family.startingGeneTree);
    for (auto &pairMapping: mappings.getMap()) {
      auto &species = pairMapping.second;
      if (speciesStringToSpeciesId.find(species) == speciesStringToSpeciesId.end()) {
        speciesStringToSpeciesId.insert({species, speciesIdToSpeciesString.size()});
        speciesIdToSpeciesString.push_back(species);
        
      }
    }
  }
  unsigned int speciesNumber = speciesIdToSpeciesString.size(); 
  std::vector<double> nullDistances(speciesNumber, 0.0);
  distanceMatrix = DistanceMatrix(speciesNumber, nullDistances);  
  DistanceMatrix distanceDenominator(speciesNumber, nullDistances);
    
  for (auto &family: families) {
    GeneSpeciesMapping mappings;
    mappings.fill(family.mappingFile, family.startingGeneTree);
    std::ifstream reader(family.startingGeneTree);
    std::string geneTreeStr;
    while (std::getline(reader, geneTreeStr)) {
      PLLUnrootedTree geneTree(geneTreeStr, false);
      geneDistancesFromGeneTree(geneTree, 
          mappings,
          speciesStringToSpeciesId,
          distanceMatrix,
          distanceDenominator,
          minMode,
          reweight,
          useBL,
          useBootstrap,
          ustar);
    }
  }
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      if (i == j) {
        distanceMatrix[i][j] = 0.0;
      }
      if (0.0 != distanceDenominator[i][j]) {
        distanceMatrix[i][j] /= distanceDenominator[i][j];
      }
    }
  }
}
