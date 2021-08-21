#include "BME.hpp"
#include <NJ/MiniNJ.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>

static void fillDistancesRec(pll_unode_t *currentNode, 
    const StringToUint &speciesStringToSpeciesId,
    double currentDistance,
    std::vector<double> &distances)
{
  currentDistance += 1.0;
  if (!currentNode->next) {
    // leaf
    std::string label(currentNode->label);
    distances[speciesStringToSpeciesId.at(label)] 
      = currentDistance;
    return;
  }
  fillDistancesRec(currentNode->next->back, 
      speciesStringToSpeciesId,
      currentDistance, 
      distances);
  fillDistancesRec(currentNode->next->next->back, 
      speciesStringToSpeciesId,
      currentDistance, 
      distances);
}

static void fillSpeciesDistances(const PLLUnrootedTree &speciesTree,
    const StringToUint &speciesStringToSpeciesId,
    DistanceMatrix &speciesDistanceMatrix)
{
  for (auto leaf: speciesTree.getLeaves()) {
    std::string label = leaf->label;
    auto id = speciesStringToSpeciesId.at(label);
    fillDistancesRec(leaf->back,
        speciesStringToSpeciesId,
        0.0,
        speciesDistanceMatrix[id]);
  }
}

double BME::computeBME(const PLLUnrootedTree &speciesTree,
    const Families &families)
{
  bool minMode = true;
  bool reweight = false;
  bool ustar = false;
  DistanceMatrix geneDistanceMatrix;
  std::vector<std::string> speciesIdToSpeciesString;
  StringToUint speciesStringToSpeciesId;
  MiniNJ::computeDistanceMatrix(families,
      minMode, 
      reweight,
      ustar,
      geneDistanceMatrix,
      speciesIdToSpeciesString,
      speciesStringToSpeciesId);

  unsigned int N = speciesIdToSpeciesString.size();
  std::vector<double> nullDistances(N, 0.0);
  DistanceMatrix speciesDistanceMatrix = DistanceMatrix(N, 
      nullDistances);  
  fillSpeciesDistances(speciesTree, 
      speciesStringToSpeciesId,
      speciesDistanceMatrix);
  double res = 0.0;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      assert(speciesDistanceMatrix[i][j] != 0.0);
      double weight = pow(2.0, speciesDistanceMatrix[i][j]);
      res += geneDistanceMatrix[i][j] / weight;
      //Logger::info << geneDistanceMatrix[i][j] << "-" << speciesDistanceMatrix[i][j] << " ";
    }
    //Logger::info << std::endl;
  }
  return res;
}

