#include "NeighborJoining.hpp"
#include <IO/Logger.hpp>

using Cherry = std::pair<unsigned int, unsigned int>;
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

static Position findMinPosition(const DistanceMatrix &distanceMatrix)
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

static Position findConstrainedMinPosition(
    const std::set<Cherry> &cherries,
    const std::vector<unsigned int> &constrainIndexToMatrixIndex,
    const DistanceMatrix &distanceMatrix
    )
{
  assert(cherries.size());
  double minDistance = std::numeric_limits<double>::max();
  const Position invalidPosition = {-1, -1};
  Position minPosition = invalidPosition;
  for (auto &cherry: cherries) {
    auto i = constrainIndexToMatrixIndex[cherry.first];
    auto j = constrainIndexToMatrixIndex[cherry.second];
    if (i >=j) {
      std::swap(i, j);
    }
    if (minDistance > distanceMatrix[i][j]) {
      minPosition.first = i;
      minPosition.second = j;
      minDistance = distanceMatrix[i][j];
    }
  }
  assert(minPosition != invalidPosition);
  return minPosition;
}

static void updateConstrainData(
    unsigned int speciesEntryToUpdate,
    unsigned int speciesEntryToRemove,
    double blLeft, 
    double blRight,
    std::set<Cherry> &constrainCherries,
    std::unordered_set<unsigned int> &constrainLeaves,
    std::vector<unsigned int> &constrainIndexToMatrixIndex,
    std::vector<unsigned int> &matrixIndexToConstrainIndex,
    PLLRootedTree &constrainTree)
{
  auto constrainIndexToUpdate = 
    matrixIndexToConstrainIndex[speciesEntryToUpdate];
  auto constrainIndexToRemove = 
    matrixIndexToConstrainIndex[speciesEntryToRemove];

  // Update branch lengthes
  constrainTree.getNode(constrainIndexToUpdate)->length = blLeft;
  constrainTree.getNode(constrainIndexToRemove)->length = blRight;

  // Erase outdated values
  auto cherriesOldSize = constrainCherries.size();
  auto leavesOldSize = constrainLeaves.size();
  Cherry toRemove(constrainIndexToUpdate, constrainIndexToRemove);
  constrainCherries.erase(toRemove);
  std::swap(toRemove.first, toRemove.second);
  constrainCherries.erase(toRemove);
  constrainLeaves.erase(constrainIndexToUpdate);
  constrainLeaves.erase(constrainIndexToRemove);
  assert(cherriesOldSize == constrainCherries.size() + 1);
  assert(leavesOldSize == constrainLeaves.size() + 2);

  // compute all indices
  auto matrixIndexToUpdate = 
    constrainIndexToMatrixIndex[constrainIndexToUpdate];
  assert(matrixIndexToUpdate == speciesEntryToUpdate);
  auto parentNode = constrainTree.getParent(constrainIndexToUpdate);
  auto parentConstrainIndex = parentNode->node_index;
  // update mappings
  constrainIndexToMatrixIndex[parentConstrainIndex] = matrixIndexToUpdate;
  matrixIndexToConstrainIndex[matrixIndexToUpdate] = parentConstrainIndex;
  // update with new values
  constrainLeaves.insert(parentConstrainIndex);
  auto neighborConstrainIndex = 
    constrainTree.getNeighbor(parentConstrainIndex)->node_index;
  if (constrainLeaves.find(neighborConstrainIndex) != 
      constrainLeaves.end()) {
    // neighborNode and parentNode form a new cherry
    constrainCherries.insert({parentConstrainIndex, neighborConstrainIndex});
  }
}


std::unique_ptr<PLLRootedTree> NeighborJoining::applyNJ(
    DistanceMatrix distanceMatrix,
    std::vector<std::string> speciesIdToSpeciesString,
    StringToUint speciesStringToSpeciesId,
    PLLRootedTree *constrainTree)
{

  /*
   * For the constrained branch length estimation:
   * - matrixIndex refers to the index in distanceMatrix
   * - constrainIndex refers to the index of a species in 
   *   the constrainTree (the species is accessible with
   *   constrainTree->getNode(constrainIndex)
   * - constrainLeaves are all current "leaves", where a leaf
   *   is a node that still has to be merged (it can correspond 
   *   to an internal node in the original tree). They are
   *   represented with constrain indices
   * - constrainCherries are all current species cherries
   *   to merge, indexed with constrain indices
   */
  std::set<Cherry> constrainCherries;
  std::unordered_set<unsigned int> constrainLeaves;
  std::vector<unsigned int> constrainIndexToMatrixIndex;
  std::vector<unsigned int> matrixIndexToConstrainIndex;
  if (constrainTree) {
    constrainIndexToMatrixIndex.resize(constrainTree->getNodeNumber());
    matrixIndexToConstrainIndex.resize(constrainTree->getNodeNumber());
    for (auto leaf: constrainTree->getLeaves()) {
      auto constrainIndex = leaf->node_index;
      auto speciesLabel = std::string(leaf->label);
      auto matrixIndex = speciesStringToSpeciesId[speciesLabel];
      constrainIndexToMatrixIndex[constrainIndex] = matrixIndex;
      matrixIndexToConstrainIndex[matrixIndex] = constrainIndex;
      constrainLeaves.insert(constrainIndex);
      auto neighborNode = constrainTree->getNeighbor(constrainIndex);
      if (!neighborNode->left) {
        // left and neighbor form a cherry
        auto neighborIndex = neighborNode->node_index;
        if (constrainLeaves.find(neighborIndex) == constrainLeaves.end()) {
          // first time we find the cherry
          constrainCherries.insert({constrainIndex, neighborIndex});
        }
      }
    }
  }

  unsigned int speciesNumber = speciesIdToSpeciesString.size();
  std::string subtree;
  // Start Neighbor Joining iterations
  for (unsigned int step = 0; step < speciesNumber - 1; ++step) {
    unsigned int dmSize = speciesNumber - step;
    auto Q = getQ(distanceMatrix, dmSize);    
    Position minPosition;
    if (constrainTree) {
      minPosition = findConstrainedMinPosition(constrainCherries,
          constrainIndexToMatrixIndex,
          Q);
    } else {
      minPosition = findMinPosition(Q);
    }
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
    bl1 = std::max(0.0, bl1);
    bl2 = std::max(0.0, bl2);
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
    
    if (constrainTree && step != speciesNumber - 2) {
      updateConstrainData(speciesEntryToUpdate,
          speciesEntryToRemove,
          bl1,
          bl2,
          constrainCherries,
          constrainLeaves,
          constrainIndexToMatrixIndex,
          matrixIndexToConstrainIndex,
          *constrainTree);
    }
  }
  std::string newick = subtree + ";";
  auto res = std::make_unique<PLLRootedTree>(newick, false); 
  
  return res;
}





