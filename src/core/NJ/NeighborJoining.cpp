#include "NeighborJoining.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLUnrootedTree.hpp>

typedef std::vector<double> Count;
typedef std::vector< std::vector<double> > DistanceMatrix;
static const double invalidDouble = std::numeric_limits<double>::infinity();

static double l1(const Count &c1, const Count &c2) {
  double res = 0.0;
  assert(c1.size() == c2.size());
  for (unsigned int i = 0; i < c1.size(); ++i) {
    res += fabs(c1[i] - c2[i]);
  }
  return res;
}

static double distance(const Count &c1, const Count &c2) {
  return l1(c1, c2);
}

static double getIfOk(double value) {
  return (value != invalidDouble) ? value : 0.0;
}

static DistanceMatrix getQ(DistanceMatrix &distanceMatrix)
{
  DistanceMatrix Q = distanceMatrix;
  unsigned int n = distanceMatrix.size();
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      double sum = 0.0;
      if (distanceMatrix[i][j] == invalidDouble) {
        sum = invalidDouble;
      } else {
        sum = (n - 2) * getIfOk(distanceMatrix[i][j]);
        for (unsigned int k = 0; k < n; ++k) {
          sum -= getIfOk(distanceMatrix[i][k]) + getIfOk(distanceMatrix[j][k]); 
        }
      }
      Q[i][j] = sum;
    }
  }
  return Q;
}

typedef std::pair<unsigned int, unsigned int> Position;

Position findMinPosition(const DistanceMatrix &distanceMatrix)
{
  double minDistance = std::numeric_limits<double>::max();
  const Position invalidPosition = {-1, -1};
  Position minPosition = invalidPosition;
  for (unsigned int i = 0; i < distanceMatrix.size(); ++i) {
    for (unsigned int j = i + 1; j < distanceMatrix.size(); ++j) {
      
      if (minDistance > distanceMatrix[i][j]) {
        minPosition.first = i;
        minPosition.second = j;
        minDistance = distanceMatrix[i][j];
      }
    }
  }
  assert(minPosition != invalidPosition);
  return minPosition;
}


static std::unique_ptr<PLLRootedTree> applyNJ(DistanceMatrix &distanceMatrix,
    std::vector<std::string> &speciesIdToSpeciesString,
    std::unordered_map<std::string, unsigned int> &speciesStringToSpeciesId)
{
  unsigned int speciesNumber = speciesIdToSpeciesString.size();
  std::string subtree;
  for (unsigned int step = 0; step < speciesNumber - 1; ++step) {
    auto Q = getQ(distanceMatrix);    
    auto minPosition = findMinPosition(Q);
    auto p1 = minPosition.first;
    auto p2 = minPosition.second;
    auto speciesEntryToUpdate = p1;
    auto speciesEntryToRemove = p2;
    // update new values:
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      double newDistance = 0.5 * (distanceMatrix[i][p1] + distanceMatrix[i][p2] - distanceMatrix[p1][p2]);
      distanceMatrix[speciesEntryToUpdate][i] = newDistance;
      distanceMatrix[i][speciesEntryToUpdate] = newDistance;
    }
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      distanceMatrix[speciesEntryToRemove][i] = invalidDouble;
      distanceMatrix[i][speciesEntryToRemove] = invalidDouble;
    }
    subtree = "(" + speciesIdToSpeciesString[minPosition.first] + ","
      + speciesIdToSpeciesString[minPosition.second] + ")";
    speciesIdToSpeciesString[speciesEntryToUpdate] = subtree;
    speciesIdToSpeciesString[speciesEntryToRemove] = std::string("NULL");
    speciesStringToSpeciesId[subtree] = speciesEntryToUpdate;
  }
  std::string newick = subtree + ";";
  return std::make_unique<PLLRootedTree>(newick, false); 
}

std::unique_ptr<PLLRootedTree> NeighborJoining::countProfileNJ(const Families &families)
{
  std::vector<Count> speciesToCountVector;
  std::vector<std::string> speciesIdToSpeciesString;
  std::unordered_map<std::string, unsigned int> speciesStringToSpeciesId;
  unsigned int familiesNumber = families.size();

  for (unsigned int i = 0; i < familiesNumber; ++i) {
    GeneSpeciesMapping mappings;
    mappings.fill(families[i].mappingFile, families[i].startingGeneTree);
    for (auto &pairMapping: mappings.getMap()) {
      auto &species = pairMapping.second;
      if (speciesStringToSpeciesId.find(species) == speciesStringToSpeciesId.end()) {
        speciesToCountVector.push_back(Count(familiesNumber, 0));
        speciesStringToSpeciesId.insert({species, speciesIdToSpeciesString.size()});
        speciesIdToSpeciesString.push_back(species);
        
      }
      speciesToCountVector[speciesStringToSpeciesId[species]][i]++;
    }
  }
  
  unsigned int speciesNumber = speciesToCountVector.size();
  std::vector<double> nullDistances(speciesNumber, 0.0);
  DistanceMatrix distanceMatrix(speciesNumber, nullDistances);
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      distanceMatrix[i][j] = distance(speciesToCountVector[i],
          speciesToCountVector[j]);
    }
  }
  return applyNJ(distanceMatrix, speciesIdToSpeciesString, speciesStringToSpeciesId);
}
 

void fillDistancesRec(pll_unode_t *currentNode, 
    bool useBL,
    double currentDistance,
    std::vector<double> &distances)
{
  if (!currentNode->next) {
    // leaf
    distances[currentNode->node_index] = currentDistance;
    return;
  }
  if (useBL) {
    currentDistance += currentNode->length;
  } else {
    currentDistance += 1.0;
  }
  fillDistancesRec(currentNode->next->back, useBL, currentDistance, distances);
  fillDistancesRec(currentNode->next->next->back, useBL, currentDistance, distances);
} 


void geneDistancesFromGeneTree(PLLUnrootedTree &geneTree,
    GeneSpeciesMapping &mapping,
    std::unordered_map<std::string, unsigned int> &speciesStringToSpeciesId,
    DistanceMatrix &distances,
    DistanceMatrix &distancesDenominator,
    bool minMode = true,
    bool normalize = false,
    bool useBL = false)
{
  unsigned int speciesNumber = distances.size();
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
    fillDistancesRec(leafNode->back, useBL, 0.0, leafDistances[leafNode->node_index]);
    fillDistancesRec(leafNode, useBL, 0.0, leafDistances[leafNode->node_index]);
  }

  // fill species distance matrices
  for (auto gene1: leaves) {
    auto gid1 = gene1->node_index;
    auto spid1 = geneIdToSpeciesId[gid1];
    for (auto gene2: leaves) {
      auto gid2 = gene2->node_index;
      auto spid2 = geneIdToSpeciesId[gid2];
      if (!minMode) {
        distances[spid1][spid2] += leafDistances[gid1][gid2];
        distancesDenominator[spid1][spid2]++;
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
  if (minMode) {
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      for (unsigned int j = 0; j < speciesNumber; ++j) {
        if (normalize) {
          distancesToAdd[i][j] /= double(leaves.size());
        }
        distances[i][j] += distancesToAdd[i][j];
        distancesDenominator[i][j] += distancesDenominatorToAdd[i][j];
      }
    }
  }
}

std::unique_ptr<PLLRootedTree> NeighborJoining::geneTreeNJ(const Families &families)
{
  std::vector<std::string> speciesIdToSpeciesString;
  std::unordered_map<std::string, unsigned int> speciesStringToSpeciesId;
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
  DistanceMatrix distanceMatrix(speciesNumber, nullDistances);
  DistanceMatrix distanceDenominator(speciesNumber, nullDistances);
  for (auto &family: families) {
    GeneSpeciesMapping mappings;
    mappings.fill(family.mappingFile, family.startingGeneTree);
    PLLUnrootedTree geneTree(family.startingGeneTree);
    geneDistancesFromGeneTree(geneTree, 
        mappings,
        speciesStringToSpeciesId,
        distanceMatrix,
        distanceDenominator);
  }
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      if (0.0 != distanceDenominator[i][j]) {
        distanceMatrix[i][j] /= distanceDenominator[i][j];
      }
      Logger::info << distanceMatrix[i][j] << " ";
    }
    Logger::info << std::endl;
  }
  
  return applyNJ(distanceMatrix, speciesIdToSpeciesString, speciesStringToSpeciesId);
}
