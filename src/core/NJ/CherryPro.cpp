#include "CherryPro.hpp"

#include <vector>
#include <memory>
#include <set>
#include "MiniNJ.hpp"
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <unordered_map>
#include <algorithm>

typedef std::vector<double> VectorDouble;
typedef std::vector<VectorDouble> MatrixDouble;
typedef std::unordered_map<std::string, int> StringToInt;
typedef std::unordered_set<int> GeneIdsSet;
typedef std::unordered_map<int, GeneIdsSet> SpeciesIdToGeneIds;
typedef std::pair<int, int> CherryProEdge;
typedef std::set<int> Clade;

static const bool CHERRY_DBG = true;
static const int PARENT_INDEX = 0;
static const int LEFT_INDEX = 2;
static const int RIGHT_INDEX = 1;
typedef std::array<int, 3> TripletInt;
typedef std::array<bool, 3> TripletBool;

static const bool USE_CHERRY_PRO_METRIC = true;
static const bool USE_WEIGHTED_CHERRY_PRO_METRIC = false;
static const bool ALWAYS_RETAG = false;

struct CherryProNode {
  bool isLeaf; 
  int sons[3];
  int geneId;
  bool speciation;
  int speciesId; // only for leaves
  bool isValid;

  CherryProNode():
    isLeaf(false),
    geneId(-1),
    speciation(true),
    speciesId(-1),
    isValid(false)
  {
  }
};

/**
 *  Unrooted gene tree for CherryPro algorithm
 */
class CherryProTree {
public:
  CherryProTree(const std::string &treeString, 
      const GeneSpeciesMapping &mapping,
      const StringToInt &speciesStrToId);
  
  std::string toNewick();

  void mergeNodesWithSameSpeciesId();
  void mergeNodesWithSpeciesId(unsigned int speciesId);
  void relabelNodesWithSpeciesId(unsigned int speciesId, 
    unsigned int newSpeciesId);

  void updateNeigborMatrix(MatrixDouble &neighborMatrix,
      MatrixDouble &denominatorMatrix);
  int coveredSpeciesNumber();
  int getCladeSize(int speciesId)  {return _speciesIdToCladeSize[speciesId];}
  std::string toString();
  void printInternalState();
  unsigned int getLeavesNumber() {
    return _leavesNumber;
  }
  std::pair<int, int> getChildrenIds(const CherryProEdge &edge);
  void findBestRootAndTag(bool removeRoot);
  int _hackIndex;
  int getRootId() {return _nodes.size() - 1;}
  int countTripletRec(int geneId, 
    const TripletInt &clade,
    TripletBool &present);
private:
  int getRootedLeftChildId(int geneId);
  int getRootedRightChildId(int geneId);
  int getRootedParentId(int geneId);
  int tagFromNode(int geneId, Clade &clade);
  
  void getAllEdgesRec(std::vector<CherryProEdge> &edges, const CherryProEdge &parentEdge);
  void getAllEdges(std::vector<CherryProEdge> &edges);
  int tagFromEdge(const CherryProEdge &edge, bool keepRoot);
  void updateDenominatorMatrixRec(int geneId, MatrixDouble &denominatorMatrix, Clade &clade);

  int getAnyValidId();
  int getAnyValidLeafId();
  int getNeighborLeaf(int nodeId);
  std::string recursiveToString(int nodeId);
  std::vector<CherryProNode> _nodes;
  SpeciesIdToGeneIds _speciesIdToGeneIds;
  std::unordered_map<int, int>  _speciesIdToCladeSize;
  unsigned int _leavesNumber;
  static int hackCounter;
};


typedef std::vector<std::shared_ptr<CherryProTree> > GeneTrees;

int CherryProTree::hackCounter = 0;

std::pair<int, int> CherryProTree::getChildrenIds(const CherryProEdge &edge)
{
  auto &node = _nodes[edge.second];
  assert(!node.isLeaf);
  if (node.sons[0] == edge.first) {
    return {node.sons[1], node.sons[2]};
  } else if (node.sons[1] == edge.first) {
    return {node.sons[0], node.sons[2]};
  } else if (node.sons[2] == edge.first) {
    return {node.sons[0], node.sons[1]};
  } else {
    assert(false);
    return {-1, -1};
  }
}


void CherryProTree::getAllEdgesRec(std::vector<CherryProEdge> &edges, const CherryProEdge &parentEdge)
{
  edges.push_back(parentEdge);
  auto &currentNode = _nodes[parentEdge.second];
  if (currentNode.isLeaf) {
    return;
  }
  auto children = getChildrenIds(parentEdge);
  getAllEdgesRec(edges, {parentEdge.second, children.first});
  getAllEdgesRec(edges, {parentEdge.second, children.second});
}

void CherryProTree::getAllEdges(std::vector<CherryProEdge> &edges)
{
  auto &anyLeaf = _nodes[getAnyValidLeafId()];
  getAllEdgesRec(edges, {anyLeaf.geneId, _nodes[anyLeaf.sons[0]].geneId});
  assert(edges.size() == 2 * _leavesNumber - 3);
}


/**
 *  Reorder the sons of node such that the son with
 *  ID parentId is set to the parent position
 */
static void setParent(CherryProNode &node, int parentId)
{
  if (node.isLeaf) {
    assert(parentId == node.sons[PARENT_INDEX]);
    return;
  }
  if (node.sons[PARENT_INDEX] == parentId) {
    // do nothing
  } else if (node.sons[(PARENT_INDEX + 1) % 3] == parentId) {
    std::swap(node.sons[(PARENT_INDEX + 1) % 3], node.sons[PARENT_INDEX]);
  } else if (node.sons[(PARENT_INDEX + 2) % 3] == parentId) {
    std::swap(node.sons[(PARENT_INDEX + 2) % 3], node.sons[PARENT_INDEX]);
  } else {
    assert(false);
  }
}


static void changeParent(CherryProNode &node, int oldParentId, int newParentId)
{
  if (node.isLeaf) {
    assert(node.sons[0] == oldParentId);
    node.sons[0] = newParentId;
    return;
  }
  for (int i = 0; i < 3; ++i) {
    if (node.sons[i] == oldParentId) {
      node.sons[i] = newParentId;
      setParent(node, newParentId);
      return;
    }
  }
  assert(false);
}

int CherryProTree::getRootedLeftChildId(int geneId)
{
  return _nodes[geneId].sons[LEFT_INDEX];
}

int CherryProTree::getRootedRightChildId(int geneId)
{
  return _nodes[geneId].sons[RIGHT_INDEX];
}

int CherryProTree::getRootedParentId(int geneId)
{
  return _nodes[geneId].sons[PARENT_INDEX];
}

static bool isIntersectionEmpty(const Clade &c1,
    const Clade &c2)
{
  auto first1 = c1.begin();
  auto first2 = c2.begin();
  auto last1 = c1.end();
  auto last2 = c2.end();
  while (first1!=last1 && first2!=last2)
  {
    if (*first1<*first2) {
      ++first1;
    } else if (*first2<*first1) {
      ++first2;
    } else {
      return false;
    }
  }
  return true;
}


/**
 *  This directly comes from astral-pro paper
 */
int CherryProTree::tagFromNode(int geneId, Clade &clade)
{
  
  auto &node = _nodes[geneId];
  if (node.isLeaf) {
    clade.insert(node.speciesId);
    return 0;
  }
  auto leftId = getRootedLeftChildId(geneId);
  auto rightId = getRootedRightChildId(geneId);
  setParent(_nodes[leftId], geneId);
  setParent(_nodes[rightId], geneId);
  Clade leftClade;
  Clade rightClade;
  auto leftScore = tagFromNode(leftId, leftClade);
  auto rightScore = tagFromNode(rightId, rightClade);
  auto score = leftScore + rightScore;
  clade.insert(leftClade.begin(), leftClade.end());
  clade.insert(rightClade.begin(), rightClade.end());
  if (isIntersectionEmpty(leftClade, rightClade)) {
    node.speciation = true;
  } else {
    node.speciation = false;
    if (leftClade == clade || rightClade == clade) {
      if (leftClade == rightClade) {
        score += 1;
      } else {
        score += 2;
      }
    } else {
      score += 3;
    }
  }
  return score;
}

int CherryProTree::tagFromEdge(const CherryProEdge &edge, bool keepRoot)
{
  // root from edge
  _nodes.push_back(CherryProNode());
  auto &root = _nodes.back();
  root.isLeaf = false;
  root.isValid = true;
  root.geneId = _nodes.size() - 1;
  root.sons[0] = -1;
  root.sons[1] = edge.first;
  root.sons[2] = edge.second;
  auto &left = _nodes[edge.first];
  auto &right = _nodes[edge.second];
  changeParent(left, right.geneId, root.geneId);
  changeParent(right, left.geneId, root.geneId);
  
  //tag from the root
  Clade clade;
  auto score = tagFromNode(root.geneId, clade);
  /*
  if (CHERRY_DBG) {
    if (keepRoot) {
      std::cout << "Score " << score << std::endl;
      std::cout << toString() << std::endl;
    }
  }
  */
  // unroot from edge
  if (!keepRoot) {
    changeParent(left, root.geneId, right.geneId);
    changeParent(right, root.geneId, left.geneId);
    _nodes.pop_back();
  }
  return score;
}

void CherryProTree::findBestRootAndTag(bool removeRoot)
{
  if (removeRoot) {
    auto &root = _nodes[getRootId()];  
    auto &left = _nodes[getRootedLeftChildId(getRootId())];
    auto &right = _nodes[getRootedRightChildId(getRootId())];
    changeParent(left, root.geneId, right.geneId);
    changeParent(right, root.geneId, left.geneId);
    _nodes.pop_back();
  }
  
  std::vector<CherryProEdge> edges;
  getAllEdges(edges);
  CherryProEdge bestEdge = {-1, -1};
  int bestScore = 10000000;
  for (auto &edge: edges) {
    int score = tagFromEdge(edge, false);
    if (score < bestScore) {
      bestEdge = edge;
      bestScore = score;
    }
  }
  assert(bestEdge.first != -1);
  assert(bestEdge.second != -1);
  tagFromEdge(bestEdge, true);
}
  
void CherryProTree::printInternalState()
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
 

static double getRatio(const MatrixDouble m, int i, int j, double bestScore)
{
  return fabs(bestScore - m[i][j]) / (bestScore + m[i][j]);
}

static bool shouldWeFight(const MatrixDouble &neighborMatrix, 
    const std::pair<int, int> &p,
    std::pair<int, int> &bestCompetingPair)
{
  return false;
  static const double RATIO_LIMIT = 0.1;
  double myScore = neighborMatrix[p.first][p.second];
  double worstRatio = 1000.0;
  for (unsigned int i = 0; i < neighborMatrix.size(); ++i) {
    auto ratio1 = getRatio(neighborMatrix, p.first, i, myScore);
    if (int(i) != p.second) {
      if (ratio1 < worstRatio) {
        worstRatio = ratio1;
        bestCompetingPair = {p.first, i};
      }
    }
    auto ratio2 = getRatio(neighborMatrix, p.second, i, myScore);
    if (int(i) != p.first) {
      if (ratio2 < worstRatio) {
        worstRatio = ratio2;
        bestCompetingPair = {p.second, i};
      }
    }
  }
  Logger::info << "Worst ratio: " << worstRatio << " from pair " << bestCompetingPair.first << " " << bestCompetingPair.second << std::endl;
  return (worstRatio < 0.1);
}


void CherryProTree::relabelNodesWithSpeciesId(unsigned int speciesId, 
    unsigned int newSpeciesId)
{
  if (_speciesIdToGeneIds.find(speciesId) == _speciesIdToGeneIds.end()) {
    return;
  }
  _speciesIdToCladeSize[newSpeciesId] += _speciesIdToCladeSize[speciesId];
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

void CherryProTree::updateDenominatorMatrixRec(int geneId, MatrixDouble &denominatorMatrix, Clade &clade)
{
  auto &node = _nodes[geneId];
  if (node.isLeaf) {
    clade.insert(node.speciesId);
    return;
  } 
  auto leftId = getRootedLeftChildId(geneId);
  auto rightId = getRootedRightChildId(geneId);
  Clade leftClade;
  Clade rightClade;
  updateDenominatorMatrixRec(leftId, denominatorMatrix, leftClade);
  updateDenominatorMatrixRec(rightId, denominatorMatrix, rightClade);
  if (node.speciation) {
    for (auto &spid1: leftClade) {
      for (auto &spid2: rightClade) {
        denominatorMatrix[spid1][spid2]++;
        denominatorMatrix[spid2][spid1]++;
      }
    }
  }
  clade.insert(leftClade.begin(), leftClade.end());
  clade.insert(rightClade.begin(), rightClade.end());
}

void CherryProTree::updateNeigborMatrix(MatrixDouble &neighborMatrixToUpdate,
      MatrixDouble &denominatorMatrixToUpdate)
{
  MatrixDouble *neighborMatrix = &neighborMatrixToUpdate;
  MatrixDouble *denominatorMatrix = &denominatorMatrixToUpdate;
  

  // First fill neighborMatrix
  for (auto &p: _speciesIdToGeneIds) {
    auto speciesId = p.first;
    const auto &geneIdSet = p.second;
    for (auto geneId: geneIdSet) {
      auto neighborGeneId = getNeighborLeaf(geneId);
      if (neighborGeneId == -1) {
        // no leaf neighbor
        continue;
      }
      auto &gene1 = _nodes[geneId];
      auto &gene2 = _nodes[neighborGeneId];
      auto spid1 = gene1.speciesId;
      auto spid2 = gene2.speciesId;
      assert(spid1 == speciesId);
      double neighborFrequency = 1.0;
      (*neighborMatrix)[spid1][spid2] += neighborFrequency;
    }
  }
  // Then fill denominatorMatrix with
  // the maximum possible number of neighbors
  // between two species
  assert(!(USE_CHERRY_PRO_METRIC && USE_WEIGHTED_CHERRY_PRO_METRIC));
  if (USE_CHERRY_PRO_METRIC) {
    Clade clade;
    updateDenominatorMatrixRec(getRootId(), *denominatorMatrix, clade);
  } else if (USE_WEIGHTED_CHERRY_PRO_METRIC) {
    Clade clade;
    int speciesNumber = denominatorMatrix->size();
    VectorDouble zeros(speciesNumber, 0.0);
    MatrixDouble maxNeighbors(speciesNumber, zeros);
    updateDenominatorMatrixRec(getRootId(), maxNeighbors, clade);
    for (auto &p: _speciesIdToGeneIds) {
      auto speciesId = p.first;
      const auto &geneIdSet = p.second;
      for (auto &p2: _speciesIdToGeneIds) {
        auto spid2 = p2.first; 
        const auto &geneIdSet2 = p2.second;
        double small = maxNeighbors[speciesId][spid2] / 2.0;
        double big = (geneIdSet.size() +  geneIdSet2.size()) / 2.0;
        (*denominatorMatrix)[speciesId][spid2] += 
          2.0  / (1.0 / big + 1.0 / small);
      }
    }
  } else {
    for (auto &p: _speciesIdToGeneIds) {
      auto speciesId = p.first;
      const auto &geneIdSet = p.second;
      for (auto &p2: _speciesIdToGeneIds) {
        auto spid2 = p2.first; 
        const auto &geneIdSet2 = p2.second;
        double small = std::min(geneIdSet.size(),  geneIdSet2.size());
        double big = std::max(geneIdSet.size(),  geneIdSet2.size());
        (*denominatorMatrix)[speciesId][spid2] += 
          2.0  / (1.0 / big + 1.0 / small);
      }
    }
  }
}
  
int CherryProTree::coveredSpeciesNumber()
{
  return _speciesIdToGeneIds.size();
}
  
void CherryProTree::mergeNodesWithSameSpeciesId()
{
  for (auto &p: _speciesIdToGeneIds) {
    mergeNodesWithSpeciesId(p.first);
  }
}


int CherryProTree::getAnyValidId()
{
  return getAnyValidLeafId();
}
 
int CherryProTree::getAnyValidLeafId()
{
  for (auto &p: _speciesIdToGeneIds) {
    for (auto id: p.second) {
      assert(_nodes[id].isLeaf);
      return id;
    }
  }
  assert(false);
  return -1;
}

 
int CherryProTree::getNeighborLeaf(int nodeId)
{
  auto &node = _nodes[nodeId];
  assert(node.isLeaf);
  assert(node.isValid);
  auto &parentNode = _nodes[node.sons[0]];
  if(parentNode.isLeaf) {
    Logger::info << "Error in " << _hackIndex << std::endl;
  }
  assert(!parentNode.isLeaf);
  for (int i = 1; i < 3; ++i) {
    auto &candidateNode = _nodes[parentNode.sons[i]];
    if (candidateNode.isLeaf && candidateNode.geneId != node.geneId) {
      return candidateNode.geneId;
    }
  }
  return -1;
}


void CherryProTree::mergeNodesWithSpeciesId(unsigned int speciesId)
{
  if (_speciesIdToGeneIds.find(speciesId) == _speciesIdToGeneIds.end()) {
    return;
  }
  auto &geneSet = _speciesIdToGeneIds[speciesId];
  GeneIdsSet geneSetCopy = geneSet;
  for (auto geneId: geneSetCopy) { 
    while (true) {
      if (getLeavesNumber() < 4) {
        // nothing to merge!
        return;
      }
      if (geneSet.find(geneId) == geneSet.end()) {
        // we already erased this gene
        break;
      }
      auto geneIdNeighbor = getNeighborLeaf(geneId);
      if (geneIdNeighbor == -1) {
        // this leaf neighbor is not a leaf, continue
        break;
      }
      if (geneSet.find(geneIdNeighbor) != geneSet.end()) {
        auto &left = _nodes[geneId];
        auto &right = _nodes[geneIdNeighbor];
        auto &parent = _nodes[left.sons[0]];
        parent.isLeaf = true;
        assert(parent.geneId != getRootId());
        // merge geneId and geneIdNeighbor
        parent.speciesId = speciesId;
        geneSet.insert(parent.geneId);
        geneSet.erase(geneId);
        geneSet.erase(geneIdNeighbor);
        left.isValid = false;
        right.isValid = false;
        geneId = parent.geneId;
        _leavesNumber--;
      } else {
        break;
      }
    }
  }
}
  

std::string CherryProTree::toString()
{
  return recursiveToString(getRootId()) + ";";
}

std::string CherryProTree::recursiveToString(int nodeId)
{
  auto &node = _nodes[nodeId];
  assert(node.isValid);
  // edge case: leaf
  if (node.isLeaf) {
    return std::to_string(node.speciesId);
  }
  auto str1 = recursiveToString(getRootedLeftChildId(nodeId));
  auto str2 = recursiveToString(getRootedRightChildId(nodeId));
  if (str1 > str2) {
    std::swap(str1, str2);
  }
  std::string res = "(";
  res += str1;
  res += ",";
  res += str2;
  res += ")";
  if (node.speciation) {
    res += "S";
  } else {
    res += "D";
  }
  return res;
}

static std::vector<int> computePLLIdToId(PLLUnrootedTree &pllTree)
{
  int currentId = 0;
  int maxPllNodeId = pllTree.getLeavesNumber()
    + pllTree.getInnerNodesNumber() * 3;
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

CherryProTree::CherryProTree(const std::string &treeString, 
      const GeneSpeciesMapping &mapping,
      const StringToInt &speciesStrToId):
  _hackIndex(hackCounter++),
  _leavesNumber(0)

{
  PLLUnrootedTree pllTree(treeString, false);
  _nodes.resize(pllTree.getLeavesNumber() * 2 - 2);
  auto pllIdToId = computePLLIdToId(pllTree);
  for (auto pllNode: pllTree.getNodes()) {
    auto geneId = pllIdToId[pllNode->node_index];
    auto &nfjNode = _nodes[geneId];
    nfjNode.sons[0] = pllIdToId[pllNode->back->node_index];
    nfjNode.geneId = geneId;
    nfjNode.isValid = true;
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
      _leavesNumber++;
    }
  }
  for (auto p: _speciesIdToGeneIds) {
    _speciesIdToCladeSize.insert({p.first, 1});
  }
  findBestRootAndTag(false);
}


static void filterGeneTrees(std::vector<std::shared_ptr<CherryProTree> > &geneTrees)
{
  auto geneTreesCopy = geneTrees;
  geneTrees.clear();
  for (auto geneTree: geneTreesCopy) {
    if (geneTree->getLeavesNumber() >= 3 && geneTree->coveredSpeciesNumber() >= 3 ) {
      geneTrees.push_back(geneTree);
    }
  }
  if (CHERRY_DBG) {
    Logger::info << "Number of gene trees after filtering: " << geneTrees.size() << std::endl;;
  }
}

int CherryProTree::countTripletRec(int geneId, 
    const TripletInt &clade,
    TripletBool &present)
{
  auto &node = _nodes[geneId];
  if (node.isLeaf) {
    for (int i = 0; i < 3; ++i) {
      present[i] = (clade[i] == node.speciesId);
    }
    return 0;
  }
  TripletBool &present1 = present;
  TripletBool present2;
  int count = 0;
  count += countTripletRec(node.sons[1], clade, present1);
  count += countTripletRec(node.sons[2], clade, present2);
  if ((present1[0] && present1[1] && present2[2])
      || (present2[0] && present2[1] && present1[2])) {
    // we found the ordered triplet!
    if (node.speciation) {
      count++;
    }
    for (int i = 0; i < 3; ++i) {
      present[i] = false;
    }
  } else if ((present1[0] && present2[1]) || (present1[1] && present2[0])) {
    if (present1[2] || present2[2]) {
      for (int i = 0; i < 3; ++i) {
        present[i] = false;
      }
    } else {
      // we found the two first elements
      present[0] = present[1] = node.speciation;
    }
  } else {
    for (int i = 0; i < 3; ++i) {
      present[i] = present1[i] || present2[i];
    }
  }
  return count;
}

static int countFight(GeneTrees &geneTrees,
    const std::pair<int, int> &p1,
    const std::pair<int, int> &p2)
{
  int third = (p1.first == p2.first || p1.second == p2.first) ? p2.second : p2.first;
  TripletInt clade = {p1.first, p1.second, third};
  int count = 0;
  TripletBool present;
  for (auto &geneTree: geneTrees) {
    count += geneTree->countTripletRec(geneTree->getRootId(), 
        clade, present);
  }
  return count;
} 

/**
 *  Return true if we decide that bestPair is
 *  really better than competingPair
 */
static bool fight(GeneTrees &geneTrees,
   const std::pair<int, int> &bestPair,
   const std::pair<int, int> &competingPair)
{
  int count1 = countFight(geneTrees, bestPair, competingPair);
  int count2 = countFight(geneTrees, competingPair, bestPair);
  Logger::info << count1 << " vs " << count2 << std::endl;
  return count1 > count2;
}

std::unique_ptr<PLLRootedTree> CherryPro::geneTreeCherryPro(const Families &families)
{
  // Init gene trees and frequency matrix
  GeneTrees geneTrees;
  StringToInt speciesStrToId;
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
      geneTrees.push_back(std::make_shared<CherryProTree>( 
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
    if (ALWAYS_RETAG)
    {
      for (auto &tree: geneTrees) {
        tree->findBestRootAndTag(true);
        tree->mergeNodesWithSameSpeciesId();
      }
      filterGeneTrees(geneTrees);
    }
    if (CHERRY_DBG) {
      for (auto &tree: geneTrees) {
        Logger::info << "Tree " << tree->_hackIndex << " " <<  tree->toString() << std::endl;
      }
    }
    
    neighborMatrix = MatrixDouble(speciesNumber, zeros);
    denominatorMatrix = MatrixDouble(speciesNumber, zeros);
    for (auto geneTree: geneTrees) {
      geneTree->updateNeigborMatrix(neighborMatrix, denominatorMatrix);
    }
    if (CHERRY_DBG) {
      Logger::info << "Neighbors: " << std::endl;
      printMatrix(neighborMatrix);
      Logger::info << "Denominators: " << std::endl;
      printMatrix(denominatorMatrix);
    }
    divideMatrix(neighborMatrix, denominatorMatrix);
    if (CHERRY_DBG) {
      Logger::info << "Frequencies: " << std::endl;
      printMatrix(neighborMatrix);
    }
    // compute the two species to join, and join them
    auto bestPairSpecies = getMaxInMatrix(neighborMatrix);
    std::pair<int, int> bestCompetingPair;
    bool doFight = shouldWeFight(neighborMatrix, bestPairSpecies, bestCompetingPair);
    if (doFight) {
      Logger::info << "Fight!!" << std::endl;
      if (!fight(geneTrees, bestPairSpecies, bestCompetingPair)) {
        Logger::info << "Lost fight, swapping!" << std::endl;
        std::swap(bestPairSpecies, bestCompetingPair);
      }
    }
    if (bestPairSpecies.first == bestPairSpecies.second) {
      // edge case when we filtered out all gene trees
      std::cout << "We filtered all gene trees, taking a random pair of species..." << std::endl;
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
      Logger::info << "Remaining gene trees: " << geneTrees.size() << std::endl;
      Logger::info << "Best pair " << bestPairSpecies.first
        << " " << bestPairSpecies.second << " with distance " << neighborMatrix[bestPairSpecies.first][bestPairSpecies.second] << std::endl;
      Logger::info << "Best pair " << speciesStr1 << " " << speciesStr2 << std::endl;
      Logger::info << speciesStr1 << std::endl << speciesStr2 << std::endl;
    }
    for (auto geneTree: geneTrees) {
      geneTree->relabelNodesWithSpeciesId(bestPairSpecies.second,
        bestPairSpecies.first);
      geneTree->mergeNodesWithSpeciesId(bestPairSpecies.first);
    }
    speciesIdToStr[bestPairSpecies.first] = std::string("(") + 
      speciesStr1 + "," + speciesStr2 + ")";
    remainingSpeciesIds.erase(bestPairSpecies.second);
  }
  std::vector<std::string> lastSpecies;
  for (auto speciesId: remainingSpeciesIds) {
    lastSpecies.push_back(speciesIdToStr[speciesId]);
  }
  assert(lastSpecies.size() == 2);
  std::string newick = "(" + lastSpecies[0] + "," + lastSpecies[1] + ");";
  if (CHERRY_DBG) {
    Logger::info << newick << std::endl;
  }
  return std::make_unique<PLLRootedTree>(newick, false); 
}



