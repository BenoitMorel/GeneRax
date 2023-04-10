#include "Cherry.hpp"

#include <vector>
#include <memory>
#include "MiniNJ.hpp"
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <unordered_map>
#include<set>
#include <util/types.hpp>

using GeneIdsSet = std::unordered_set<int>;
using SpeciesIdToGeneIds = std::unordered_map<int, GeneIdsSet>;
using Clade = std::set<int>;


static const bool CHERRY_DBG = false;
static const bool WITHOUT_CHERRY_MERGING = false;
static const bool MININJ = false;
static const int TREES_TO_DISPLAY = 0;

struct CherryNode {
  bool isLeaf; 
  int sons[3];
  int geneId;
  int speciesId; // only for leaves
  bool isValid;
  int height;
  double length;
  Clade speciesClade; // indixed with the initial
                              // species indices
};

/**
 *  Unrooted gene tree for Cherry algorithm
 */
class CherryTree {
public:
  CherryTree(const std::string &treeString, 
      const GeneSpeciesMapping &mapping,
      const StringToUint &speciesStrToId);
  
  std::string toNewick();

  void mergeNodesWithSameSpeciesId();
  void mergeNodesWithSpeciesId(unsigned int speciesId);
  void mergeNodesWithSpeciesId(unsigned int speciesId1, 
      unsigned int speciesId2,
      double &sumBL1,
      double &denBL1,
      double &sumBL2,
      double &denBL2
      );
  void relabelNodesWithSpeciesId(unsigned int speciesId, 
    unsigned int newSpeciesId);

  void updateNeigborMatrix(MatrixDouble &neighborMatrix,
      MatrixDouble &denominatorMatrix);
  int coveredSpeciesNumber();
  std::string toString();
  void printInternalState();
  unsigned int getLeafNumber() {
    return _leavesNumber;
  }

  std::pair<int, int> getChildren(int geneId, int parentId);
  void updateMiniNJMatrix(MatrixDouble &distanceMatrixToUpdate,
      MatrixDouble &denominatorMatrixToUpdate);
  int _hackIndex;
private:
  // return child geneId or -1 if not exists
  int getChildIdWithSpecies(int geneId, int speciesId);
  int mergeNeighbors(int geneId1, int geneId2);
  int getAnyValidId();
  int getNeighborLeaf(int nodeId);
  std::string recursiveToString(int nodeId, int sonToSkipId);
  void erase(int nodeId);
  void updateInternodeDistance(int geneId, 
    int parentGeneId, 
    int sourceSpeciesId,
    double currentDistance,
    MatrixDouble &speciesDistance,
    MatrixDouble &denominator,
    MatrixDouble &count);
  int getThirdIndex(CherryNode &node,
    int id1, 
    int id2);
  void connect(CherryNode &node,
    int previousId,
    int newId);

  std::vector<CherryNode> _nodes;
  SpeciesIdToGeneIds _speciesIdToGeneIds;
  unsigned int _leavesNumber;
  static int hackCounter;
};

int CherryTree::hackCounter = 0;
  
void CherryTree::printInternalState()
{
  for (auto &node: _nodes) {
    if (node.isValid) {
      Logger::info << "gid=" << node.geneId;
      if (node.isLeaf) {
        Logger::info << " spid=" << node.speciesId << "parent=" << node.sons[0] << std::endl;
      } else {
        Logger::info << " " << node.sons[0] << " "  << node.sons[1] << " " << node.sons[2] << std::endl;
      }
    }
  }
  Logger::info << "Covered species " << coveredSpeciesNumber() << std::endl;
}

static void printMatrix(const MatrixDouble m) 
{
  for (auto &v: m) {
    for (auto e: v) {
      Logger::info << e << " ";
    }
    Logger::info << std::endl;
  }

}

static void divideMatrix(MatrixDouble &m, const MatrixDouble &denom)
{
  assert(m.size() == denom.size());
  assert(m[0].size() == denom[0].size());
  for (unsigned int i = 0; i < m.size(); ++i) {
    for (unsigned int j = 0; j < m.size(); ++j) {
      if (denom[i][j] != 0.0) {
        m[i][j] /= denom[i][j];
      }
    }
  }
}



static std::pair<int, int> getMaxInMatrix(MatrixDouble &m)
{
  assert(m.size());
  std::pair<int, int> minPair = {0, 0};
  for (unsigned int i = 0; i < m.size(); ++i) {
    for (unsigned int j = 0; j < m.size(); ++j) {
      if (m[minPair.first][minPair.second] < m[i][j]) {
        minPair = {i, j};    
      }
    }
  }
  return minPair;
}
static std::pair<int, int> getMinInMatrix(MatrixDouble &m)
{
  assert(m.size());
  std::pair<int, int> minPair = {0, 0};
  for (unsigned int i = 0; i < m.size(); ++i) {
    for (unsigned int j = 0; j < m.size(); ++j) { 
      auto minValue = m[minPair.first][minPair.second];
      if (minValue == 0.0 || (minValue >= m[i][j] && m[i][j] != 0.0)) {
        minPair = {i, j};    
      }
    }
  }
  return minPair;
}
  
void CherryTree::relabelNodesWithSpeciesId(unsigned int speciesId, 
    unsigned int newSpeciesId
    )
{
  for (auto geneId: _speciesIdToGeneIds[speciesId]) {
    _nodes[geneId].height++;
  }
  for (auto geneId: _speciesIdToGeneIds[newSpeciesId]) {
    _nodes[geneId].height++;
  }
  if (_speciesIdToGeneIds.find(speciesId) == _speciesIdToGeneIds.end()) {
    return;
  }
  if (_speciesIdToGeneIds.find(newSpeciesId) == _speciesIdToGeneIds.end()) {
    _speciesIdToGeneIds.insert({newSpeciesId, GeneIdsSet()});
  }
  auto &newGeneIdSet = _speciesIdToGeneIds[newSpeciesId];
  for (auto geneId: _speciesIdToGeneIds[speciesId]) {
    auto &node = _nodes[geneId];
    node.speciesId = newSpeciesId;
    newGeneIdSet.insert(geneId);
  }
  _speciesIdToGeneIds.erase(speciesId);
}



void CherryTree::updateInternodeDistance(int geneId, 
    int parentGeneId, 
    int sourceSpeciesId,
    double currentDistance,
    MatrixDouble &speciesDistances,
    MatrixDouble &denominator,
    MatrixDouble &count)
{
  auto &node = _nodes[geneId];
  if (node.isLeaf) {
    currentDistance += node.height;
    auto &v = speciesDistances[sourceSpeciesId][node.speciesId];
    auto &d = denominator[sourceSpeciesId][node.speciesId];
    auto &c = count[sourceSpeciesId][node.speciesId];
    if (v == 0.0) {
      v = currentDistance; //std::min(3.0, currentDistance);
      d = 1.0;
      c = 1.0;
    } else {
      /*
      if (currentDistance < v) {
        v = currentDistance;
        c = 1.0;
      } else if (currentDistance == v) {
        c +=  1.0;
      }
      */
      v = std::min(currentDistance, v);
    }
    return;
  }
  auto children = getChildren(geneId, parentGeneId);
  updateInternodeDistance(children.first, geneId, sourceSpeciesId, currentDistance + 1.0, speciesDistances, denominator, count);
  updateInternodeDistance(children.second, geneId, sourceSpeciesId, currentDistance + 1.0, speciesDistances, denominator, count);
}

void CherryTree::updateMiniNJMatrix(MatrixDouble &distanceMatrixToUpdate,
      MatrixDouble &denominatorMatrixToUpdate)
{
  int speciesNumber = distanceMatrixToUpdate.size(); 
  double zero = 0.0;
  VectorDouble zeros(speciesNumber, zero);
  MatrixDouble speciesDistance(speciesNumber, zeros);
  MatrixDouble denominator(speciesNumber, zeros);
  MatrixDouble count(speciesNumber, zeros);
  for (auto &p: _speciesIdToGeneIds) {
    auto speciesId = p.first;
    const auto &geneIdSet = p.second;
    for (auto geneId: geneIdSet) {
      auto &node = _nodes[geneId];
      updateInternodeDistance(node.sons[0], geneId, speciesId, node.height, speciesDistance, denominator, count);
    }
  }
  for (int i = 0; i < speciesNumber; ++i) {
    for (int j = 0; j < speciesNumber; ++j) {
      if (i != j && speciesDistance[i][j] != zero) { 
        distanceMatrixToUpdate[i][j] += count[i][j] * speciesDistance[i][j];
        denominatorMatrixToUpdate[i][j] += count[i][j] * denominator[i][j];
      }
    }
  }  
}
  
std::pair<int, int> CherryTree::getChildren(int geneId, int parentId)
{
  auto &node = _nodes[geneId];
  if (node.sons[0] == parentId) {
    return {node.sons[1], node.sons[2]};
  } else if (node.sons[1] == parentId) {
    return {node.sons[0], node.sons[2]};
  } else if (node.sons[2] == parentId) {
    return {node.sons[0], node.sons[1]};
  } else {
    assert(false);
    return {-1, -1};
  }
}



void CherryTree::updateNeigborMatrix(MatrixDouble &neighborMatrixToUpdate,
      MatrixDouble &denominatorMatrixToUpdate)
{
  MatrixDouble *neighborMatrix = &neighborMatrixToUpdate;
  MatrixDouble *denominatorMatrix = &denominatorMatrixToUpdate;
  
  
  /*
  // settings
  const bool perFamilyWeight = false; 
  
  // init from settings
  auto speciesNumber = neighborMatrixToUpdate.size();
  bool intermediateMatrices = perFamilyWeight;
  bool deleteMatrices = false;
  if (intermediateMatrices) {
    VectorDouble zeros(speciesNumber);
    neighborMatrix = new MatrixDouble(speciesNumber, zeros);
    denominatorMatrix = new MatrixDouble(speciesNumber, zeros);
    deleteMatrices = true;
  }
  */

  for (auto &p: _speciesIdToGeneIds) {
    auto speciesId = p.first;
    // First fill neighborMatrix
    const auto &geneIdSet = p.second;
    for (auto geneId: geneIdSet) {
      auto neighborGeneId = getNeighborLeaf(geneId);
      if (neighborGeneId == -1) {
        // no leaf neighbor
        continue;
      }
      auto spid1 = _nodes[geneId].speciesId;
      auto spid2 = _nodes[neighborGeneId].speciesId;
      assert(spid1 == speciesId);
      (*neighborMatrix)[spid1][spid2] += 1.0;
    }
    // Then fill denominatorMatrix with
    // the maximum possible number of neighbors
    // between two species
    for (auto &p2: _speciesIdToGeneIds) {
      auto spid2 = p2.first; 
      const auto &geneIdSet2 = p2.second;
      double small = std::min(geneIdSet.size(),  geneIdSet2.size());
      double big = std::max(geneIdSet.size(),  geneIdSet2.size());
      (*denominatorMatrix)[speciesId][spid2] += 
        2.0  / (1.0 / big + 1.0 / small);
        
      //std::min(geneIdSet.size(),  geneIdSet2.size());
    }
  }
 
  /*
  // update according to settings
  if (perFamilyWeight) {
    for (unsigned int i = 0; i < speciesNumber; ++i) {
      for (unsigned int j = 0; j < speciesNumber; ++j) {
        if ((*denominatorMatrix)[i][j] != 0.0) {
          neighborMatrixToUpdate[i][j] += (*neighborMatrix)[i][j] / (*denominatorMatrix)[i][j];
          neighborMatrixToUpdate[i][j] += 1.0;
        }
      }
    }
  }

  if (deleteMatrices) {
    delete denominatorMatrix;
    delete neighborMatrix;
  }  
  */
}
  
int CherryTree::coveredSpeciesNumber()
{
  return _speciesIdToGeneIds.size();
}
  
void CherryTree::mergeNodesWithSameSpeciesId()
{
  for (auto &p: _speciesIdToGeneIds) {
    mergeNodesWithSpeciesId(p.first);
  }
}


int CherryTree::getAnyValidId()
{
  for (auto &p: _speciesIdToGeneIds) {
    for (auto id: p.second) {
      return id;
    }
  }
  assert(false);
  return -1;
}
 
int CherryTree::getNeighborLeaf(int nodeId)
{
  auto &node = _nodes[nodeId];
  assert(node.isLeaf);
  assert(node.isValid);
  auto &parentNode = _nodes[node.sons[0]];
  if(parentNode.isLeaf) {
    Logger::info << "Error in " << _hackIndex << std::endl;
  }
  assert(!parentNode.isLeaf);
  for (int i = 0; i < 3; ++i) {
    auto &candidateNode = _nodes[parentNode.sons[i]];
    if (candidateNode.isLeaf && candidateNode.geneId != node.geneId) {
      return candidateNode.geneId;
    }
  }
  return -1;
}


int CherryTree::mergeNeighbors(int geneId1, int geneId2)
{
  auto &node1 = _nodes[geneId1];
  auto &node2 = _nodes[geneId2];
  auto &parent = _nodes[node1.sons[0]];
  assert(node1.isLeaf && node2.isLeaf);
  assert(node1.sons[0] == node2.sons[0]);
  //assert(node1.speciesId == node2.speciesId);
  for (int i = 0; i < 3; ++i) {
    if (parent.sons[i] != geneId1 &&
        parent.sons[i] != geneId2) {
      parent.sons[0] = parent.sons[i];
      break;
    }
  }
  parent.height = 0;
  parent.isLeaf = true;
  parent.speciesId = node1.speciesId;
  auto &geneSet = _speciesIdToGeneIds[node1.speciesId];
  geneSet.insert(parent.geneId);
  _leavesNumber++;
  parent.speciesClade.insert(node1.speciesClade.begin(), node1.speciesClade.end());
  parent.speciesClade.insert(node2.speciesClade.begin(), node2.speciesClade.end());
  erase(node1.geneId);
  erase(node2.geneId);
  return parent.geneId;
}

int CherryTree::getChildIdWithSpecies(int geneId, int speciesId)
{
  auto &node = _nodes[geneId];
  if (node.isLeaf) {
    return -1;
  }
  for (int i = 0; i < 3; ++i) {
    auto &child = _nodes[node.sons[i]];
    if (child.isLeaf && child.speciesId == speciesId) {
      return child.geneId;
    }
  }
  return -1;
} 

int CherryTree::getThirdIndex(CherryNode &node,
    int id1, 
    int id2)
{
  for (int i = 0; i < 3; ++i) {
    if (node.sons[i] != id1 && node.sons[i] != id2) {
      return i;
    }
  }
  assert(false);
  return -1;
}

void CherryTree::connect(CherryNode &node,
    int previousId,
    int newId)
{
  for (int i = 0; i < 3; ++i) {
    if (node.sons[i] == previousId) {
      node.sons[i] = newId;
      return;
    }
  }
  assert(false);
}


void CherryTree::erase(int nodeId)
{
  auto &node = _nodes[nodeId];
  node.isValid = false;
  if (node.isLeaf) {
    _leavesNumber--;
    auto &geneSet = _speciesIdToGeneIds[node.speciesId];
    geneSet.erase(node.geneId);
  }
  node.speciesClade.clear();
}

void CherryTree::mergeNodesWithSpeciesId(unsigned int speciesId)
{
  if (_speciesIdToGeneIds.find(speciesId) == _speciesIdToGeneIds.end()) {
    return;
  }
  auto &geneSet = _speciesIdToGeneIds[speciesId];
  GeneIdsSet geneSetCopy = geneSet;
  // MERGE NEIGHBORS
  for (auto geneId: geneSetCopy) { 
    while (true) {
      if (getLeafNumber() < 4) {
        // nothing to merge!
        return;
      }
      if (geneSet.find(geneId) == geneSet.end()) {
        // we already erased this gene
        break;
      }
      auto geneIdNeighbor = getNeighborLeaf(geneId);
      if (geneIdNeighbor != -1) { // CHERRY CASE
        if (_nodes[geneIdNeighbor].speciesId == int(speciesId)) {
          // merge
          geneId = mergeNeighbors(geneId, geneIdNeighbor);
        } else {
          break;
        }
      } else {
        break;
      }
    }
  }
}
  
void CherryTree::mergeNodesWithSpeciesId(unsigned int speciesId1, 
      unsigned int speciesId2,
      double &sumBL1,
      double &denBL1,
      double &sumBL2,
      double &denBL2)
{
  if (_speciesIdToGeneIds.find(speciesId1) == _speciesIdToGeneIds.end()) {
    return;
  }
  if (_speciesIdToGeneIds.find(speciesId2) == _speciesIdToGeneIds.end()) {
    return;
  }
  auto &geneSet1 = _speciesIdToGeneIds[speciesId1];
  auto &geneSet2 = _speciesIdToGeneIds[speciesId2];
  GeneIdsSet genesToMerge = geneSet1;
  genesToMerge.insert(geneSet2.begin(), geneSet2.end());
  // MERGE NEIGHBORS
  for (auto geneId: genesToMerge) { 
    while (true) {
      if (getLeafNumber() < 4) {
        // nothing to merge!
        return;
      }
      if (geneSet1.find(geneId) == geneSet1.end()
          && geneSet2.find(geneId) == geneSet2.end()) {
        // we already erased this gene
        break;
      }
      auto geneIdNeighbor = getNeighborLeaf(geneId);
      
      if (geneIdNeighbor != -1) { // CHERRY CASE
        auto neighborSpid = _nodes[geneIdNeighbor].speciesId;
        auto spid = _nodes[geneId].speciesId;
        if (neighborSpid == int(speciesId1)) {
          if (spid != neighborSpid) {
            sumBL1 += _nodes[geneIdNeighbor].length;
            sumBL2 += _nodes[geneId].length;
            denBL1 += 1.0;
            denBL2 += 1.0;
          }
          geneId = mergeNeighbors(geneId, geneIdNeighbor);
        } else if (neighborSpid == int(speciesId2)) {
          if (spid != neighborSpid) {
            sumBL1 += _nodes[geneId].length;
            sumBL2 += _nodes[geneIdNeighbor].length;
            denBL1 += 1.0;
            denBL2 += 1.0;
          }
          geneId = mergeNeighbors(geneId, geneIdNeighbor);
        } else {
          break;
        }
      } else {
        break;
      }
    }
  }
}
  

std::string CherryTree::toString()
{
  return recursiveToString(getAnyValidId(), -1) + ";";
}

std::string CherryTree::recursiveToString(int nodeId, int sonToSkipId)
{
  auto &node = _nodes[nodeId];
  assert(node.isValid);
  // edge case: do not start the recursion from a leaf
  if (sonToSkipId == -1 && node.isLeaf) {
    return recursiveToString(node.sons[0], -1);
  }
  // edge case: leaf
  if (node.isLeaf) {
    auto res = std::to_string(node.speciesId);
    res += "{";
    for (auto &e: node.speciesClade) {
      res += std::to_string(e);
      res += ",";
    }
    res.back() = '}';
    return res;
  }
  std::string res = "(";
  bool firstNodeWritten = false;
  for (int i = 0; i < 3; ++i) {
    if (sonToSkipId != node.sons[i]) {
      res += recursiveToString(node.sons[i], node.geneId);
      if ((!firstNodeWritten || sonToSkipId == -1) && i != 2) {
        res += ",";
        firstNodeWritten = true;
      }
    }
  }
  res += ")";
  return res;
}

static std::vector<int> computePLLIdToId(PLLUnrootedTree &pllTree)
{
  int currentId = 0;
  int maxPllNodeId = pllTree.getLeafNumber()
    + pllTree.getInnerNodeNumber() * 3;
  std::vector<int> pllIdToId(maxPllNodeId, -1);
  for (auto pllNode: pllTree.getNodes()) {
    if (pllIdToId[pllNode->node_index] == -1) {
      pllIdToId[pllNode->node_index] = currentId;
      if (pllNode->next) {
        pllIdToId[pllNode->next->node_index] = currentId;
        pllIdToId[pllNode->next->next->node_index] = currentId;
      }
      currentId++;
    }
  }
  return pllIdToId;
}

CherryTree::CherryTree(const std::string &treeString, 
      const GeneSpeciesMapping &mapping,
      const StringToUint &speciesStrToId):
  _hackIndex(hackCounter++),
  _leavesNumber(0)
{
  PLLUnrootedTree pllTree(treeString, false);
  _nodes.resize(pllTree.getLeafNumber() * 2 - 2);
  auto pllIdToId = computePLLIdToId(pllTree);
  for (auto pllNode: pllTree.getNodes()) {
    auto geneId = pllIdToId[pllNode->node_index];
    auto &nfjNode = _nodes[geneId];
    nfjNode.sons[0] = pllIdToId[pllNode->back->node_index];
    nfjNode.geneId = geneId;
    nfjNode.isValid = true;
    nfjNode.height = 0;
    nfjNode.length = pllNode->length;
    if (pllNode->next) {
      // internal node
      nfjNode.isLeaf = false;
      nfjNode.sons[1] = pllIdToId[pllNode->next->back->node_index];
      nfjNode.sons[2] = pllIdToId[pllNode->next->next->back->node_index];
      nfjNode.speciesId = -1;
    } else {
      // leaf node
      nfjNode.isLeaf = true;
      auto species = mapping.getSpecies(pllNode->label);
      nfjNode.speciesId = speciesStrToId.at(species);
      if (_speciesIdToGeneIds.find(nfjNode.speciesId) == 
          _speciesIdToGeneIds.end()) {
        _speciesIdToGeneIds.insert({nfjNode.speciesId, GeneIdsSet()});
      }
      _speciesIdToGeneIds[nfjNode.speciesId].insert(nfjNode.geneId);
      nfjNode.speciesClade.insert(nfjNode.speciesId);
      _leavesNumber++;
    }
  }
}


static void filterGeneTrees(std::vector<std::shared_ptr<CherryTree> > &geneTrees)
{
  auto geneTreesCopy = geneTrees;
  geneTrees.clear();
  for (auto geneTree: geneTreesCopy) {
    if (geneTree->getLeafNumber() >= 4 && geneTree->coveredSpeciesNumber() > 2 ) {
      geneTrees.push_back(geneTree);
    }
  }
  if (CHERRY_DBG) {
    Logger::info << "Number of gene trees after filtering: " << geneTrees.size() << std::endl;;
  }
}

std::unique_ptr<PLLRootedTree> Cherry::geneTreeCherry(const Families &families)
{
  // Init gene trees and frequency matrix
  std::vector<std::shared_ptr<CherryTree> > geneTrees;
  StringToUint speciesStrToId;
  std::vector<std::string> speciesIdToStr;
  
  // fill the structure that map speciesStr <-> speciesId
  // and create gene trees mapped with the species IDs
  for (auto &family: families) {
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    auto coveredSpecies = mapping.getCoveredSpecies();
    for (auto &species: mapping.getCoveredSpecies()) {
      if (speciesStrToId.find(species) == speciesStrToId.end()) {
        speciesStrToId.insert({species, speciesIdToStr.size()});
        speciesIdToStr.push_back(species);
      }
    }
    std::ifstream reader(family.startingGeneTree);
    std::string line;
    while (std::getline(reader, line)) {
      geneTrees.push_back(std::make_shared<CherryTree>( 
          line, mapping, speciesStrToId));
    }
  }
  Logger::info << "Loaded " << geneTrees.size() << " gene trees" << std::endl;
  unsigned int speciesNumber = speciesStrToId.size();
  std::unordered_set<int> remainingSpeciesIds;
  for (unsigned int i = 0; i < speciesNumber; ++i) {
    remainingSpeciesIds.insert(i);
  }
  filterGeneTrees(geneTrees);
  for (auto geneTree: geneTrees) {
    geneTree->mergeNodesWithSameSpeciesId();
  }
  // main loop of the algorithm
  VectorDouble zeros(speciesNumber, 0.0);
  MatrixDouble neighborMatrix(speciesNumber, zeros);
  MatrixDouble denominatorMatrix(speciesNumber, zeros);
  for (unsigned int i = 0; i < speciesNumber - 2; ++i) {
    if (CHERRY_DBG) {
      Logger::info << std::endl;
      Logger::info << "*******************************" << std::endl;
      Logger::info << "Remaining species: " << remainingSpeciesIds.size() << std::endl;
      Logger::info << "Species mappings:" << std::endl;
      for (auto spid: remainingSpeciesIds) {
        Logger::info << "  " << spid << "\t" << speciesIdToStr[spid] << std::endl;
      }
    }
      // filter out gene trees that do not hold information
    filterGeneTrees(geneTrees);
    // merge identical tips in cherries
    // compute the frequency matrix
    if (CHERRY_DBG) {
      for (int i = 0; i < std::min(int(geneTrees.size()), TREES_TO_DISPLAY); ++i) {
        Logger::info << "Tree " << geneTrees[i]->_hackIndex << " " <<  geneTrees[i]->toString() << std::endl;
      }
    }
    
    neighborMatrix = MatrixDouble(speciesNumber, zeros);
    denominatorMatrix = MatrixDouble(speciesNumber, zeros);
    for (auto geneTree: geneTrees) {
      if (MININJ) {
        geneTree->updateMiniNJMatrix(neighborMatrix, denominatorMatrix);
      } else {
        geneTree->updateNeigborMatrix(neighborMatrix, denominatorMatrix);
      }
    }
    /*
    if (CHERRY_DBG) {
      Logger::info << "Neighbors: " << std::endl;
      printMatrix(neighborMatrix);
      Logger::info << "Denominators: " << std::endl;
      printMatrix(denominatorMatrix);
    }
    */
    divideMatrix(neighborMatrix, denominatorMatrix);
    if (CHERRY_DBG) {
      Logger::info << "Frequencies: " << std::endl;
      printMatrix(neighborMatrix);
    }
    if (WITHOUT_CHERRY_MERGING) {
      if (!MININJ) {
        for (auto &c: neighborMatrix) {
          for (auto &d: c) {
            d = -d;
          }
        }
      }
      return MiniNJ::applyNJ(neighborMatrix, speciesIdToStr, speciesStrToId);
    }
    // compute the two species to join, and join them
    std::pair<int, int> bestPairSpecies = MININJ ? 
      getMinInMatrix(neighborMatrix) : getMaxInMatrix(neighborMatrix);
    if (bestPairSpecies.first == bestPairSpecies.second) {
      // edge case when we filtered out all gene trees
      assert(geneTrees.size() == 0);
      bestPairSpecies = {-1, -1};
      for (auto speciesId: remainingSpeciesIds) {
        if (bestPairSpecies.first == -1) {
          bestPairSpecies.first = speciesId;
          continue;
        }
        if (bestPairSpecies.second == -1) {
          bestPairSpecies.second = speciesId;
          break;
        }
      }
      assert(bestPairSpecies.second != -1);
    }
    std::string speciesStr1 = speciesIdToStr[bestPairSpecies.first];
    std::string speciesStr2 = speciesIdToStr[bestPairSpecies.second];
    if (CHERRY_DBG) {
      Logger::info << "Best pair " << bestPairSpecies.first
        << " " << bestPairSpecies.second << " with distance " << neighborMatrix[bestPairSpecies.first][bestPairSpecies.second] << std::endl;
      Logger::info << "Best pair " << speciesStr1 << " " << speciesStr2 << std::endl;
    }
    double sumBL1 = 0.0;
    double denBL1 = 0.0;
    double sumBL2 = 0.0;
    double denBL2 = 0.0;
    for (auto geneTree: geneTrees) {
      geneTree->mergeNodesWithSpeciesId(bestPairSpecies.first,
          bestPairSpecies.second,
          sumBL1,
          denBL1,
          sumBL2,
          denBL2);
      geneTree->relabelNodesWithSpeciesId(bestPairSpecies.second,
        bestPairSpecies.first);
      //geneTree->mergeNodesWithSpeciesId(bestPairSpecies.first);
    }
    double bl1 = sumBL1 / denBL1;
    double bl2 = sumBL2 / denBL2;
    if (bl1 != bl1) {
      bl1 = 0.0;
    }
    if (bl2 != bl2) {
      bl2 = 0.0;
    }
    speciesIdToStr[bestPairSpecies.first] = std::string("(") + 
      speciesStr1 + ":" + std::to_string(bl1) + "," + 
      speciesStr2 + ":" + std::to_string(bl2) + ")";
    remainingSpeciesIds.erase(bestPairSpecies.second);
  }
  std::vector<std::string> lastSpecies;
  std::vector<unsigned int> lastSpeciesId;
  for (auto speciesId: remainingSpeciesIds) {
    lastSpecies.push_back(speciesIdToStr[speciesId]);
    lastSpeciesId.push_back(speciesId);
  }
  assert(lastSpecies.size() == 2);
  double sumBL1 = 0.0;
  double denBL1 = 0.0;
  double sumBL2 = 0.0;
  double denBL2 = 0.0;
  for (auto geneTree: geneTrees) {
    geneTree->mergeNodesWithSpeciesId(lastSpeciesId[0],
        lastSpeciesId[1],
        sumBL1,
        denBL1,
        sumBL2,
        denBL2);
  }
  double bl1 = sumBL1 / denBL1;
  double bl2 = sumBL2 / denBL2;
  if (bl1 != bl1) {
    bl1 = 0.0;
  }
  if (bl2 != bl2) {
    bl2 = 0.0;
  }
  std::string newick = "(" + lastSpecies[0] + ":" + std::to_string(bl1) +
                       "," + lastSpecies[1] + ":" + std::to_string(bl2) + ");";
  Logger::info << newick << std::endl;
  return std::make_unique<PLLRootedTree>(newick, false); 
}



