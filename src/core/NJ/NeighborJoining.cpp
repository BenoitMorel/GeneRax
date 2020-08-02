#include "NeighborJoining.hpp"
#include <IO/Logger.hpp>

static const double invalidDouble = std::numeric_limits<double>::infinity();

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

std::unique_ptr<PLLRootedTree> NeighborJoining::applyNJ(
    DistanceMatrix &distanceMatrix,
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





