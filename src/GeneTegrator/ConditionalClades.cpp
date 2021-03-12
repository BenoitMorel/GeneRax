  
#include "ConditionalClades.hpp"
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>

#include <fstream>
#include <iostream>

struct TreeWraper {
  std::shared_ptr<PLLUnrootedTree> tree;
  bool operator ==(const TreeWraper &other) const
  {
    return PLLUnrootedTree::areIsomorphic(*tree, *(other.tree));
  }
};

namespace std
{
  template<>
    struct hash<TreeWraper>
    {
      size_t
        operator()(const TreeWraper& tree) const
        {
          return tree.tree->getUnrootedTreeHash();
        }
    };
}

using WeightedTrees = std::unordered_map<TreeWraper, unsigned int>;



static CCPClade getComplementary(const CCPClade &clade,
    const CCPClade &subclade)
{
  CCPClade complementary(subclade);
  complementary.negate();
  complementary &= clade;
  return complementary;
}

static void addClade(CID cid,
    CladeCounts &cladeCounts,
    unsigned int count)
{
  auto cladeCountsIt = cladeCounts.find(cid);
  if (cladeCountsIt == cladeCounts.end()) {
    cladeCounts.insert({cid, count});
  } else {
    cladeCountsIt->second += count;
  }
}

static void addSubclade(CID cid,
    CID subcladeCID,
    SubcladeCounts &subcladeCounts,
    unsigned int count)
{
  addClade(subcladeCID, subcladeCounts[cid], count);
}
   
void printClade(const CCPClade &clade, 
    const std::vector<std::string> &idToLeaf)
{
  std::cerr << "{";
  bool first = true;
  for (unsigned int i = 0; i < clade.size(); ++i) {
    if (clade[i]) {
      if (!first) {
        std::cerr << ",";
      }
      std::cerr << idToLeaf[i];
      first = false;
    }
  }
  std::cerr << "}";
}





void readTrees(const std::string &newickFile,
    WeightedTrees &weightedTrees,
    unsigned int &inputTrees,
    unsigned int &uniqueInputTrees)
{
  std::ifstream infile(newickFile);
  std::string line;
  while (std::getline(infile, line)) {
    TreeWraper wraper;
    wraper.tree = std::make_shared<PLLUnrootedTree>(line, false);
    auto it = weightedTrees.find(wraper);
    if (it == weightedTrees.end()) {
      weightedTrees.insert({wraper, 1});
    } else {
      it->second++;
    }
    inputTrees++;
  }
  uniqueInputTrees = weightedTrees.size();
}


static void firstPass(const WeightedTrees &weightedTrees,
    const std::unordered_map<std::string, unsigned int> &leafToId,
    CladeToCID &cladeToCID,
    CIDToClade &cidToClade,
    CIDToLeaf &cidToLeaf,
    std::vector<std::vector<pll_unode_t*> > &postOrderNodes
    )
{
  auto &anyTree = *(weightedTrees.begin()->first.tree);
  CCPClade emptyClade(anyTree.getLeavesNumber(), false);
  CCPClade fullClade(anyTree.getLeavesNumber(), true);
  auto leafNumber = anyTree.getLeavesNumber();
  std::unordered_set<CCPClade> unorderedClades;
  unorderedClades.insert(fullClade);
  std::vector<pll_unode_t *> nodes;
  std::vector<char> buffer;
  for (auto pair: weightedTrees) {
    auto &tree = *(pair.first.tree);
    std::vector<CCPClade> nodeIndexToClade(tree.getDirectedNodesNumber(), 
        emptyClade);
    postOrderNodes.emplace_back(tree.getPostOrderNodes());
    for (auto node: postOrderNodes.back()) {
      auto nodeIndex = node->node_index;
      auto &clade = nodeIndexToClade[nodeIndex];
      if (!node->next) { // leaf
        auto id = leafToId.at(node->label);
        clade.set(id);
      } else {
        auto &leftClade = nodeIndexToClade[node->next->back->node_index];
        auto &rightClade = nodeIndexToClade[node->next->next->back->node_index];
        clade = leftClade | rightClade;
      }
      unorderedClades.insert(clade);
    }
  }
  OrderedClades orderedClades;
  for (auto &clade: unorderedClades) {
    orderedClades.insert(clade);
  }
  cidToClade.clear();
  cladeToCID.clear();
  for (auto it = orderedClades.begin(); it != orderedClades.end(); ++it) {
    auto &clade = *it;
    unsigned int CID = cidToClade.size();
    cidToClade.push_back(clade);
    cladeToCID[clade] = CID;
  }
  for (auto pair: leafToId) {
    CCPClade clade(leafNumber, false);
    clade.set(pair.second);
    cidToLeaf[cladeToCID[clade]] = pair.first;
  }
}



ConditionalClades::ConditionalClades(const std::string &newickFile):
  _skip(false),
  _inputTrees(0),
  _uniqueInputTrees(0)
{

  std::unordered_map<std::string, unsigned int> leafToId;
  CCPClade emptyClade;
  CCPClade fullClade;
  WeightedTrees weightedTrees;
  readTrees(newickFile, weightedTrees, _inputTrees, _uniqueInputTrees); 
  if (_inputTrees == _uniqueInputTrees) {
    _skip = true;
    return;
  }
  auto &anyTree = *(weightedTrees.begin()->first.tree);
  for (auto leaf: anyTree.getLabels()) {
    leafToId.insert({leaf, leafToId.size()});
    _idToLeaf.push_back(leaf);
  }
  emptyClade = CCPClade(anyTree.getLeavesNumber(), false);
  fullClade = CCPClade(anyTree.getLeavesNumber(), true);
  
  std::vector<std::vector<pll_unode_t*> > postOrderNodes;
  firstPass(weightedTrees, 
      leafToId,
      _cladeToCID,
      _CIDToClade,
      _CIDToLeaf,
      postOrderNodes
      );
  CladeCounts cladeCounts;
  SubcladeCounts subcladeCounts(_CIDToClade.size());
  auto fullCladeCID = _cladeToCID[fullClade];
  std::vector<CCPClade> nodeIndexToClade(anyTree.getDirectedNodesNumber(), emptyClade); 
  std::vector<CID> nodeIndexToCID(anyTree.getDirectedNodesNumber());
  // root clade
  addClade(fullCladeCID, cladeCounts, 
      _inputTrees * anyTree.getDirectedNodesNumber());
  unsigned int weightedTreeIndex = 0;
  for (auto pair: weightedTrees) {
    auto &tree = *(pair.first.tree);
    auto treeCount = pair.second;
    for (auto node: postOrderNodes[weightedTreeIndex]) {
      auto nodeIndex = node->node_index;
      auto &clade = nodeIndexToClade[nodeIndex];
      CID cid;
      if (!node->next) { // leaf
        auto id = leafToId[node->label];
        clade = emptyClade;
        clade.set(id);
        cid = _cladeToCID[clade];
      } else {
        auto &leftClade = nodeIndexToClade[node->next->back->node_index];
        auto &rightClade = nodeIndexToClade[node->next->next->back->node_index];
        auto leftCID = nodeIndexToCID[node->next->back->node_index];
        auto rightCID = nodeIndexToCID[node->next->next->back->node_index];
        clade = leftClade | rightClade;
        cid = _cladeToCID[clade];
        if (leftCID < rightCID) {
          auto leftCID = _cladeToCID[leftClade];
          addSubclade(cid, leftCID, subcladeCounts, treeCount);
        } else {
          auto rightCID = _cladeToCID[rightClade];
          addSubclade(cid, rightCID, subcladeCounts, treeCount);
        }
      }
      nodeIndexToCID[node->node_index] = cid;
      addClade(cid, cladeCounts, treeCount);
      
      // todo should we check that a clade/complementary is only
      // added once to the root??
      addSubclade(fullCladeCID, cid, subcladeCounts, treeCount);
    }
    weightedTreeIndex++;
  }
  _fillCCP(cladeCounts, subcladeCounts);
}
  
void ConditionalClades::_fillCCP(CladeCounts &cladeCounts,
      SubcladeCounts &subcladeCounts)
{
  _allCladeSplits.clear();
  auto cladesNumber = _CIDToClade.size();
  _allCladeSplits.resize(cladesNumber);
  for (unsigned int cid = 0; cid < cladesNumber; ++cid) {
    auto &clade = _CIDToClade[cid];
    auto cladeCount = cladeCounts[cid];
    auto &cladeSplits = _allCladeSplits[cid];
    if (subcladeCounts[cid].size()) {
      // internal clade
      double sumFrequencies = 0.0;
      for (auto &subcladeCount: subcladeCounts[cid]) {
        auto CIDLeft = subcladeCount.first;
        auto cladeLeft = _CIDToClade[CIDLeft];
        auto cladeRight = getComplementary(clade, cladeLeft);
        auto CIDRight = _cladeToCID[cladeRight];
        double frequency = double(subcladeCount.second) / double(cladeCount); 
        sumFrequencies += frequency;
        CladeSplit split;
        split.parent = cid;
        split.left = CIDLeft;
        split.right = CIDRight;
        split.frequency = frequency;
        cladeSplits.push_back(split);
      }
      // Check that frequencies sum to one
      assert(fabs(1.0 - sumFrequencies) < 0.000001);
    }
  }
}

void ConditionalClades::printContent() const
{

  for (unsigned int CID = 0; CID < _allCladeSplits.size(); ++CID) {
    auto &clade = _CIDToClade[CID];
    auto &cladeSplits = _allCladeSplits[CID];
    printClade(clade, _idToLeaf);
    std::cout << std::endl;
    double sumFrequencies = 0.0;
    for (auto &cladeSplit: cladeSplits) {
      std::cerr << "  ";
      printClade(_CIDToClade[cladeSplit.left], _idToLeaf);
      std::cerr << "|";
      printClade(_CIDToClade[cladeSplit.right], _idToLeaf);
      sumFrequencies += cladeSplit.frequency;
      std::cerr << " -> " << cladeSplit.frequency << std::endl;
    }
  }
}
  
bool ConditionalClades::isLeaf(CID cid) const
{
  return getCladeSplits(cid).size() == 0;
}
  
std::string ConditionalClades::getLeafLabel(CID cid) const
{
  assert(isLeaf(cid));
  return _CIDToLeaf.at(cid);
}
  
unsigned int ConditionalClades::getRootsNumber() const
{
  return _allCladeSplits.back().size() / 2; 
}

