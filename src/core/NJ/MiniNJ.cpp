#include "MiniNJ.hpp"




static double getIfOk(double value) {
  return (value != invalidDouble) ? value : 0.0;
}

void printDistanceMatrix(const DistanceMatrix &d) {
  for (auto &l: d) {
    for (auto elem: l) {
      Logger::info << elem << "\t";
    }
    Logger::info << std::endl;
  }
  Logger::info << std::endl;
}

static DistanceMatrix getQ(DistanceMatrix &distanceMatrix, unsigned int dmSize)
{
  DistanceMatrix Q = distanceMatrix;
  unsigned int n = dmSize;
  for (unsigned int i = 0; i < distanceMatrix.size(); ++i) {
    for (unsigned int j = 0; j < distanceMatrix.size(); ++j) {
      double sum = 0.0;
      if (distanceMatrix[i][j] == invalidDouble) {
        sum = invalidDouble;
      } else {
        sum = (n - 2) * getIfOk(distanceMatrix[i][j]);
        for (unsigned int k = 0; k < distanceMatrix.size(); ++k) {
          sum -= getIfOk(distanceMatrix[i][k]) + getIfOk(distanceMatrix[j][k]); 
        }
      }
      Q[i][j] = sum;
    }
  }
  return Q;
}

using Position = std::pair<unsigned int, unsigned int>;

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


std::unique_ptr<PLLRootedTree> MiniNJ::applyNJ(DistanceMatrix &distanceMatrix,
    std::vector<std::string> &speciesIdToSpeciesString,
    StringToUint &speciesStringToSpeciesId)
{
  unsigned int speciesNumber = speciesIdToSpeciesString.size();
  std::string subtree;
  for (unsigned int step = 0; step < speciesNumber - 1; ++step) {
    unsigned int dmSize = speciesNumber - step;
    auto Q = getQ(distanceMatrix, dmSize);    
    auto minPosition = findMinPosition(Q);
    auto p1 = minPosition.first;
    auto p2 = minPosition.second;
    auto speciesEntryToUpdate = p1;
    auto speciesEntryToRemove = p2;
    // update new values:
    DistanceMatrix copy = distanceMatrix;
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      double newDistance = 0.5 * (copy[i][p1] + copy[i][p2] - copy[p1][p2]);
      distanceMatrix[speciesEntryToUpdate][i] = newDistance;
      distanceMatrix[i][speciesEntryToUpdate] = newDistance;
    }
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      distanceMatrix[speciesEntryToRemove][i] = invalidDouble;
      distanceMatrix[i][speciesEntryToRemove] = invalidDouble;
    }
    double bl1 = 0.0;
    for (unsigned int k = 0; k < speciesNumber; ++k) {
      if (copy[p1][k] != invalidDouble) {
        bl1 += copy[p1][k] - copy[p2][k];
      }
    }
    if (dmSize > 2) {
      bl1 /= double(dmSize - 2);
    }
    bl1 += copy[p1][p2];
    bl1 *= 0.5;
    double bl2 = copy[p1][p2] - bl1;

    Logger::info << bl1 << " " << bl2 << std::endl;

    subtree = "(" + 
      speciesIdToSpeciesString[minPosition.first] 
      + ":" + std::to_string(bl1)
      + ","
      + speciesIdToSpeciesString[minPosition.second] 
      + ":" + std::to_string(bl2)
      + ")";
    speciesIdToSpeciesString[speciesEntryToUpdate] = subtree;
    speciesIdToSpeciesString[speciesEntryToRemove] = std::string("NULL");
    speciesStringToSpeciesId[subtree] = speciesEntryToUpdate;
  }
  std::string newick = subtree + ";";
  return std::make_unique<PLLRootedTree>(newick, false); 
}

void fillDistancesRec(pll_unode_t *currentNode, 
    bool useBL,
    bool useBootstrap,
    double currentDistance,
    std::vector<double> &distances)
{
  if (useBL) {
    currentDistance += currentNode->length;
  } else if (useBootstrap) {
    double bootstrapValue = (nullptr == currentNode->label) ? 0.0 : std::atof(currentNode->label);
    currentDistance += bootstrapValue;
  } else {
    if (currentNode->length > 0.0000011) {
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


void geneDistancesFromGeneTree(PLLUnrootedTree &geneTree,
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
    //fillDistancesRec(leafNode, useBL, useBootstrap, 0.0, leafDistances[leafNode->node_index]);
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


std::unique_ptr<PLLRootedTree> MiniNJ::geneTreeNJ(const Families &families, bool minAlgo, bool ustarAlgo, bool reweight)
{
  std::vector<std::string> speciesIdToSpeciesString;
  StringToUint speciesStringToSpeciesId;
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
    
  bool minMode = minAlgo;
  bool useBL = false;
  bool useBootstrap = false;
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
          ustarAlgo);
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
   // Logger::info << std::endl;
  }
  return applyNJ(distanceMatrix, speciesIdToSpeciesString, speciesStringToSpeciesId);
}
