  
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





void readTrees(const std::string &inputFile,
    WeightedTrees &weightedTrees,
    unsigned int &inputTrees,
    unsigned int &uniqueInputTrees,
    int maxSamples)
{
  std::ifstream infile(inputFile);
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(infile, line)) {
    lines.push_back(line);
  }
  if (maxSamples >= 0 && static_cast<int>(lines.size()) > maxSamples) {
    // only get a subset of the trees
    auto copy = lines;
    lines.clear();
    for (int i = 0; i < maxSamples; ++i) {
      lines.push_back(copy[i * (copy.size() / maxSamples)]);
    }
  }
  for (auto line: lines) {
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

// extract all clades and map them to clade IDs, which 
// are sorted (child clades appear before the parent clades)
// for each weighted tree, also stores its node in post order
static void extractClades(const WeightedTrees &weightedTrees,
    const std::unordered_map<std::string, unsigned int> &leafToId,
    CladeToCID &cladeToCID,
    CIDToClade &cidToClade,
    CIDToLeaf &cidToLeaf,
    std::vector<std::vector<corax_unode_t*> > &postOrderNodes
    )
{
  auto &anyTree = *(weightedTrees.begin()->first.tree);
  CCPClade emptyClade(anyTree.getLeavesNumber(), false);
  CCPClade fullClade(anyTree.getLeavesNumber(), true);
  auto leafNumber = anyTree.getLeavesNumber();
  std::vector<corax_unode_t *> nodes;
  std::vector<char> buffer;
  // iterate over all trees of the tree distribution
  OrderedClades orderedClades;
  orderedClades.insert(fullClade);
  for (auto pair: weightedTrees) {
    auto &tree = *(pair.first.tree);
    std::vector<CCPClade> nodeIndexToClade(tree.getDirectedNodesNumber(), 
        emptyClade);
    postOrderNodes.emplace_back(tree.getPostOrderNodes());
    // traverse all directed nodes of the current tree, and compute
    // the corresponding clade. Store them into nodeIndexToClade to 
    // reuse the child subclades in the post order traversal
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
      orderedClades.insert(clade);
    }
  }
  cidToClade.clear();
  cladeToCID.clear();
  // fill cidToClade and cladeToCID
  // The order is important because we want to make sure that a parent
  // clade always comes after its child clades in the CID ordering 
  for (auto it = orderedClades.begin(); it != orderedClades.end(); ++it) {
    auto &clade = *it;
    unsigned int CID = cidToClade.size();
    cidToClade.push_back(clade);
    cladeToCID[clade] = CID;
  }
  // so far we have only added the non-trivial clades. Now we add
  // the trivial clades 
  for (auto pair: leafToId) {
    CCPClade clade(leafNumber, false);
    clade.set(pair.second);
    cidToLeaf[cladeToCID[clade]] = pair.first;
  }
}


// fill cladeCounts (number of times each clade occurs) and
// subcladeCOunts (for each clade, the number of times each subclade
// occures)
static void fillCladeCounts(const WeightedTrees &weightedTrees,
    unsigned int inputTrees,
    const CladeToCID &cladeToCID,
    const std::unordered_map<std::string, unsigned int> &leafToId,
    const std::vector<std::vector<corax_unode_t*> > &postOrderNodes,
    CladeCounts &cladeCounts,
    SubcladeCounts &subcladeCounts)
{
  auto &anyTree = *(weightedTrees.begin()->first.tree);
  auto emptyClade = CCPClade(anyTree.getLeavesNumber(), false);
  auto fullClade = CCPClade(anyTree.getLeavesNumber(), true);
  auto fullCladeCID = cladeToCID.at(fullClade);
  std::vector<CCPClade> nodeIndexToClade(anyTree.getDirectedNodesNumber(), emptyClade); 
  std::vector<CID> nodeIndexToCID(anyTree.getDirectedNodesNumber());
  // root clade
  addClade(fullCladeCID, cladeCounts, 
      inputTrees * anyTree.getDirectedNodesNumber());
  unsigned int weightedTreeIndex = 0;
  for (auto pair: weightedTrees) {
    auto treeCount = pair.second;
    for (auto node: postOrderNodes.at(weightedTreeIndex)) {
      auto nodeIndex = node->node_index;
      auto &clade = nodeIndexToClade.at(nodeIndex);
      CID cid;
      if (!node->next) { // leaf
        auto id = leafToId.at(node->label);
        clade = emptyClade;
        clade.set(id);
        cid = cladeToCID.at(clade);
      } else {
        auto &leftClade = nodeIndexToClade[node->next->back->node_index];
        auto &rightClade = nodeIndexToClade[node->next->next->back->node_index];
        auto leftCID = nodeIndexToCID[node->next->back->node_index];
        auto rightCID = nodeIndexToCID[node->next->next->back->node_index];
        clade = leftClade | rightClade;
        cid = cladeToCID.at(clade);
        if (leftCID < rightCID) {
          auto leftCID = cladeToCID.at(leftClade);
          addSubclade(cid, leftCID, subcladeCounts, treeCount);
        } else {
          auto rightCID = cladeToCID.at(rightClade);
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

}



ConditionalClades::ConditionalClades(const std::string &inputFile,
    bool fromBinary, 
    int maxSamples):
  _inputTrees(0),
  _uniqueInputTrees(0)
{
  if (fromBinary) {
    unserialize(inputFile);
    return;
  }
  std::unordered_map<std::string, unsigned int> leafToId;
  CCPClade emptyClade;
  CCPClade fullClade;
  WeightedTrees weightedTrees;
  readTrees(inputFile, weightedTrees, _inputTrees, _uniqueInputTrees, maxSamples); 
  auto &anyTree = *(weightedTrees.begin()->first.tree);
  for (auto leaf: anyTree.getLabels()) {
    leafToId.insert({leaf, leafToId.size()});
    _idToLeaf.push_back(leaf);
  }
  emptyClade = CCPClade(anyTree.getLeavesNumber(), false);
  fullClade = CCPClade(anyTree.getLeavesNumber(), true);
  
  std::vector<std::vector<corax_unode_t*> > postOrderNodes;
  // first pass to get all clades and assign them a CID
  extractClades(weightedTrees, 
      leafToId,
      _cladeToCID,
      _CIDToClade,
      _CIDToLeaf,
      postOrderNodes
      );


  CladeCounts cladeCounts;
  SubcladeCounts subcladeCounts(_CIDToClade.size());
  // second pass to count the number of occurence of each clade,
  // and for each clade, the number of occurences of its subclades
  fillCladeCounts(weightedTrees,
      _inputTrees,
      _cladeToCID,
      leafToId,
      postOrderNodes,
      cladeCounts,
      subcladeCounts);
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

static void serializeUInt(unsigned int v,
    std::ostream &os)
{
  os.write(reinterpret_cast<char *>(&v),sizeof(unsigned int));
}

static unsigned int unserializeUInt(std::istream &is)
{
  unsigned int res;
  is.read(reinterpret_cast<char *>(&res),sizeof(unsigned int));
  return res;
}

static void serializeDouble( double v,
    std::ostream &os)
{
  os.write(reinterpret_cast<char *>(&v),sizeof(double));
}

static double unserializeDouble(std::istream &is)
{
  double res;
  is.read(reinterpret_cast<char *>(&res),sizeof(double));
  return res;
}

static void serializeString(const std::string &str,
    std::ostream &os)
{
  unsigned int size = str.size();
  serializeUInt(str.size(), os);
  os.write(str.c_str(), size*sizeof(char) );
}

static std::string unserializeString(std::istream &is)
{
  unsigned int size = unserializeUInt(is);
  std::vector<char> temp(size);
  is.read(reinterpret_cast<char *>(&temp[0]), size * sizeof(char));
  return std::string (temp.begin(),temp.end());
}

static void serializeCCPClade(const CCPClade &clade,
    std::ostream &os)
{
  unsigned int size = clade.size();
  serializeUInt(size, os);
  auto &data = clade.getInternalBuffer();
  unsigned int dataSize = data.size();
  serializeUInt(dataSize, os);
  os.write(reinterpret_cast<const char *>(&data[0]),
      data.size() * sizeof(genesis::utils::Bitvector::IntType));
}

static CCPClade unserializeCCPClade(std::istream &is)
{
  auto size = unserializeUInt(is);
  CCPClade res(size);
  auto &data = res.getInternalBuffer();
  auto dataSize = unserializeUInt(is);
  data.resize(dataSize);
  is.read(reinterpret_cast<char *>(&data[0]), 
      data.size() * sizeof(genesis::utils::Bitvector::IntType));
  return res;
}

static void serializeCladeSplit(const CladeSplit &split,
    std::ostream &os)
{
  serializeUInt(split.parent, os); 
  serializeUInt(split.left, os); 
  serializeUInt(split.right, os); 
  serializeDouble(split.frequency, os); 
}

static CladeSplit unserializeCladeSplit(std::ifstream &is)
{
  CladeSplit res;
  res.parent = unserializeUInt(is);
  res.left = unserializeUInt(is);
  res.right = unserializeUInt(is);
  res.frequency = unserializeDouble(is);
  return res;
}

void ConditionalClades::serialize(const std::string &outputFile)
{
  std::ofstream os(outputFile, std::ios::binary);
  serializeUInt(_inputTrees, os);
  serializeUInt(_uniqueInputTrees, os);
  // _idToLeaf
  serializeUInt(_idToLeaf.size(), os);
  for (const auto &str: _idToLeaf) {
    serializeString(str, os);
  }
  // _CIDToLeaf
  serializeUInt(_CIDToLeaf.size(), os);
  for (auto it: _CIDToLeaf) {
    serializeUInt(it.first, os);
    serializeString(it.second, os);
  }
  // _cladeToCID
  serializeUInt(_cladeToCID.size(), os);
  for (auto it: _cladeToCID) {
    serializeCCPClade(it.first, os);
    serializeUInt(it.second, os);
  }
  // _CIDToClade
  serializeUInt(_CIDToClade.size(), os);
  for (auto &clade: _CIDToClade) {
    serializeCCPClade(clade, os);
  }
  // _allCladeSplits
  serializeUInt(_allCladeSplits.size(), os);
  for (auto &cladeSplits: _allCladeSplits) {
    serializeUInt(cladeSplits.size(), os);
    for (auto &cladeSplit: cladeSplits) {
      serializeCladeSplit(cladeSplit, os); 
    }
  }

  Logger::info << _cladeToCID.size() << std::endl;
}

void ConditionalClades::unserialize(const std::string &inputFile)
{
  std::ifstream is(inputFile, std::ios::binary);
  _inputTrees = unserializeUInt(is);
  _uniqueInputTrees = unserializeUInt(is);
  // _idToLeaf
  _idToLeaf.resize(unserializeUInt(is));
  for (unsigned int i = 0; i < _idToLeaf.size(); ++i) {
    _idToLeaf[i] = unserializeString(is);
  }
  // _CIDToLeaf
  auto cidToLeafSize = unserializeUInt(is);
  for (unsigned int i = 0; i < cidToLeafSize; ++i) {
    auto cid = unserializeUInt(is);
    auto leaf = unserializeString(is);
    _CIDToLeaf.insert({cid, leaf});
  }
  // _cladeToCID
  auto cladeToCIDSize = unserializeUInt(is);
  for (unsigned int i = 0; i < cladeToCIDSize; ++i) {
    auto clade = unserializeCCPClade(is);
    auto CID = unserializeUInt(is);
    _cladeToCID.insert({clade, CID});
  }
  // _CIDToClade
  auto cidToCladeSize = unserializeUInt(is);
  _CIDToClade.resize(cidToCladeSize);
  for (unsigned int i = 0; i < cidToCladeSize; ++i) {
    _CIDToClade[i] = unserializeCCPClade(is);
  }
  // _allCladeSplits
  auto allCladesSplitSize = unserializeUInt(is);
  _allCladeSplits.resize(allCladesSplitSize);
  for (unsigned int i = 0; i < allCladesSplitSize; ++i) {
    auto &cladeSplits = _allCladeSplits[i];
    auto cladeSplitsSize = unserializeUInt(is);
    cladeSplits.resize(cladeSplitsSize);
    for (unsigned int j = 0; j < cladeSplitsSize; ++j) {
      cladeSplits[j] = unserializeCladeSplit(is);
    }
  }
  Logger::info << _cladeToCID.size() << std::endl;
  Logger::info << std::endl;
}

void ConditionalClades::printStats() const
{
  std::cerr << "Leaves " << _CIDToLeaf.size() << " clades: " << _cladeToCID.size() <<  " unique trees " << _uniqueInputTrees << std::endl;

}
