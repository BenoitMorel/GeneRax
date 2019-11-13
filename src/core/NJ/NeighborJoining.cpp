#include "NeighborJoining.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>


typedef std::vector<double> Count;
typedef std::vector< std::vector<double> > DistanceMatrix;
static const double invalidDouble = std::numeric_limits<double>::infinity();

// for debugging
static void printDistanceMatrix(const DistanceMatrix &distanceMatrix) {
  unsigned int speciesNumber = distanceMatrix.size();
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      Logger::info << distanceMatrix[i][j] << "\t";
    }
    Logger::info << std::endl;
  }
}

// for debugging
static void printCount(const Count &count) {
  Logger::info << "[";
  for (auto c: count) {
    Logger::info << c << ", ";
  }
  
  Logger::info << "]" << std::endl;
}

static double l1(const Count &c1, const Count &c2) {
  double res = 0.0;
  assert(c1.size() == c2.size());
  for (unsigned int i = 0; i < c1.size(); ++i) {
    res += fabs(c1[i] - c2[i]);
  }
  return res;
}

static double l2(const Count &c1, const Count &c2) {
  double res = 0.0;
  assert(c1.size() == c2.size());
  for (unsigned int i = 0; i < c1.size(); ++i) {
    res += (c1[i] - c2[i]) * (c1[i] - c2[i]);
  }
  return sqrt(res);
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
    //Logger::info << "species " << speciesIdToSpeciesString[i] << ": ";
    //printCount(speciesToCountVector[i]);
    for (unsigned int j = 0; j < speciesNumber; ++j) {
      //Logger::info << "species " << speciesIdToSpeciesString[j] << ": ";
      //printCount(speciesToCountVector[i]);
      distanceMatrix[i][j] = distance(speciesToCountVector[i],
          speciesToCountVector[j]);
    }
  }

  std::string subtree;
  for (unsigned int step = 0; step < speciesNumber - 1; ++step) {
    auto Q = getQ(distanceMatrix);    
    /*
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      Logger::info << speciesIdToSpeciesString[i] << "\t";
    }
    Logger::info << std::endl;
    Logger::info << "DistanceMatrix: " << std::endl;
    printDistanceMatrix(distanceMatrix);
    Logger::info << "Q:" << std::endl;
    printDistanceMatrix(Q);
    */
    auto minPosition = findMinPosition(Q);
    auto minDistance = Q[minPosition.first][minPosition.second];
    //Logger::info << "Min position: " << minPosition.first << " " << minPosition.second << " " << minDistance << std::endl;
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
    //Logger::info << subtree << std::endl;
    speciesIdToSpeciesString[speciesEntryToUpdate] = subtree;
    speciesIdToSpeciesString[speciesEntryToRemove] = std::string("NULL");
    speciesStringToSpeciesId[subtree] = speciesEntryToUpdate;
    
  }
  std::string newick = subtree + ";";
  return std::make_unique<PLLRootedTree>(newick, false); 
}

