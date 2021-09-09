#include "MiniBME.hpp"
#include <NJ/MiniNJ.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>


static DistanceMatrix getNullMatrix(unsigned int N)
{
  std::vector<double> nullDistances(N, 0.0);
  return DistanceMatrix(N, nullDistances); 
}


static void fillDistancesRec(pll_unode_t *currentNode, 
    const StringToUint &speciesStringToSpeciesId,
    double currentDistance,
    std::vector<double> &distances,
    std::vector<bool> *belongsToPruned)
{
  if (!belongsToPruned || (*belongsToPruned)[currentNode->node_index]) {
    currentDistance += 1.0;
  }
  if (!currentNode->next) {
    // leaf
    if (!belongsToPruned || (*belongsToPruned)[currentNode->node_index]) {
      std::string label(currentNode->label);
      distances[speciesStringToSpeciesId.at(label)] 
        = currentDistance;
    }
    return;
  }
  fillDistancesRec(currentNode->next->back, 
      speciesStringToSpeciesId,
      currentDistance, 
      distances,
      belongsToPruned);
  fillDistancesRec(currentNode->next->next->back, 
      speciesStringToSpeciesId,
      currentDistance, 
      distances,
      belongsToPruned);
}

static void fillSpeciesDistances(const PLLUnrootedTree &speciesTree,
    const StringToUint &speciesStringToSpeciesId,
    DistanceMatrix &speciesDistanceMatrix,
    std::vector<bool> *belongsToPruned = nullptr)
{
  for (auto leaf: speciesTree.getLeaves()) {
    if (belongsToPruned && !(*belongsToPruned)[leaf->node_index]) {
      continue;
    }
    std::string label = leaf->label;
    auto id = speciesStringToSpeciesId.at(label);
    fillDistancesRec(leaf->back,
        speciesStringToSpeciesId,
        0.0,
        speciesDistanceMatrix[id],
        belongsToPruned);
  }
}

static DistanceMatrix getPrunedSpeciesMatrix(const PLLUnrootedTree &speciesTree,
  const StringToUint &speciesStringToSpeciesId,
  const std::unordered_set<std::string> &coverage)
{
  unsigned int N = speciesTree.getLeavesNumber();
  DistanceMatrix res = getNullMatrix(N);
  std::vector<bool> hasChildren(speciesTree.getDirectedNodesNumber(), false);
  std::vector<bool> belongsToPruned(speciesTree.getDirectedNodesNumber(), false);
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto index = node->node_index;
    if (node->next) {
      auto leftIndex = node->next->back->node_index;
      auto rightIndex = node->next->next->back->node_index;
      hasChildren[index] = hasChildren[leftIndex] || hasChildren[rightIndex];
      belongsToPruned[index] = hasChildren[leftIndex] && hasChildren[rightIndex];
    } else {
      bool isCovered = coverage.find(std::string(node->label)) != coverage.end();
      hasChildren[index] = isCovered;
      belongsToPruned[index] = isCovered;
    }
  }
  fillSpeciesDistances(speciesTree,
      speciesStringToSpeciesId,
      res,
      &belongsToPruned);
  return res;
}
    




MiniBME::MiniBME(const PLLUnrootedTree &speciesTree,
    const Families &families,
    bool pruneMode):
  _prune(pruneMode)
{
  bool minMode = true;
  bool reweight = false;
  bool ustar = false;
  if (_prune) {
    // fill species-spid mappings
    for (auto leaf: speciesTree.getLeaves()) {
      std::string label(leaf->label);
      _speciesStringToSpeciesId.insert({label, _speciesIdToSpeciesString.size()});
      _speciesIdToSpeciesString.push_back(label);
    }
    PerCoreGeneTrees::getPerCoreFamilies(families, _perCoreFamilies);
    _geneDistanceMatrices.resize(_perCoreFamilies.size());
    _geneDistanceDenominators.resize(_perCoreFamilies.size());
    unsigned int speciesNumber = speciesTree.getLeavesNumber();
    std::vector<double> nullDistances(speciesNumber, 0.0);
    _perFamilyCoverage.resize(_perCoreFamilies.size());
    // map each species string to a species ID
    for (unsigned int i = 0; i < _perCoreFamilies.size(); ++i) {
      auto &family = _perCoreFamilies[i];
      _geneDistanceMatrices[i] = DistanceMatrix(speciesNumber, nullDistances);  
      _geneDistanceDenominators[i] = DistanceMatrix(speciesNumber, nullDistances);  
      GeneSpeciesMapping mappings;
      mappings.fill(family.mappingFile, family.startingGeneTree);
      for (auto species: mappings.getCoveredSpecies()) {
        _perFamilyCoverage[i].insert(species);
      }
      std::ifstream reader(family.startingGeneTree);
      std::string geneTreeStr;
      while (std::getline(reader, geneTreeStr)) {
        PLLUnrootedTree geneTree(geneTreeStr, false);
        MiniNJ::geneDistancesFromGeneTree(geneTree, 
            mappings,
            _speciesStringToSpeciesId,
            _geneDistanceMatrices[i],
            _geneDistanceDenominators[i],
            minMode,
            reweight,
            ustar);
      }
    }
  } else {  
    _geneDistanceMatrices.resize(1);
    MiniNJ::computeDistanceMatrix(families,
      minMode, 
      reweight,
      ustar,
      _geneDistanceMatrices[0],
      _speciesIdToSpeciesString,
      _speciesStringToSpeciesId);
  }
}


double MiniBME::computeBME(const PLLUnrootedTree &speciesTree)
{
  if (_prune) {
    return _computeBMEPrune(speciesTree);
  }
  unsigned int N = _speciesIdToSpeciesString.size();
  std::vector<double> nullDistances(N, 0.0);
  DistanceMatrix speciesDistanceMatrix = DistanceMatrix(N, 
      nullDistances);  
  fillSpeciesDistances(speciesTree, 
      _speciesStringToSpeciesId,
      speciesDistanceMatrix);
  double res = 0.0;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      double weight = pow(2.0, speciesDistanceMatrix[i][j]);
      res += _geneDistanceMatrices[0][i][j] / weight;
    }
  }
  return res;
}


double MiniBME::_computeBMEPrune(const PLLUnrootedTree &speciesTree)
{
  unsigned int N = _speciesIdToSpeciesString.size();
  DistanceMatrix speciesDistanceMatrix = getNullMatrix(N);
  fillSpeciesDistances(speciesTree, 
      _speciesStringToSpeciesId,
      speciesDistanceMatrix);
  double res = 0.0;
  std::vector<DistanceMatrix> prunedSpeciesMatrices(_perCoreFamilies.size(), getNullMatrix(N));
  for (unsigned int k = 0; k < _perCoreFamilies.size(); ++k) {
    auto pruned = 
    prunedSpeciesMatrices[k] = getPrunedSpeciesMatrix(speciesTree, 
          _speciesStringToSpeciesId,
          _perFamilyCoverage[k]);    
  }
  /*
  for (unsigned int k = 0; k < _perCoreFamilies.size(); ++k) {
    for (unsigned int i = 0; i < N; ++i) {
      for (unsigned int j = 0; j < N; ++j) {
        Logger::info << prunedSpeciesMatrices[k][i][j] << " ";
      }
      Logger::info << std::endl;
    }
    Logger::info << std::endl;
  } 
  */
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      double distance = 0.0;
      double denominator = 0.0;
      for (unsigned k = 0; k < _geneDistanceMatrices.size(); ++k) {
        if (0.0 == _geneDistanceDenominators[k][i][j]) {
          continue;
        }
        /*
        double correction = (speciesDistanceMatrix[i][j] / prunedSpeciesMatrices[k][i][j]);
        distance += _geneDistanceMatrices[k][i][j] * correction;
        */
        double correction = speciesDistanceMatrix[i][j] - prunedSpeciesMatrices[k][i][j];
        distance += _geneDistanceMatrices[k][i][j] + correction;
        denominator += _geneDistanceDenominators[k][i][j];
        //Logger::info << i << "-" << j << "-" << k << " " << prunedSpeciesMatrices[k][i][j] << " " <<  _geneDistanceMatrices[k][i][j] << " " << _geneDistanceDenominators[k][i][j] << std::endl;
      }
      ParallelContext::sumDouble(distance);
      ParallelContext::sumDouble(denominator);
      if (denominator != 0.0) {
        denominator *= pow(2.0, speciesDistanceMatrix[i][j]);
        res += distance / denominator;
      }
    }
  }
  return res;
}


