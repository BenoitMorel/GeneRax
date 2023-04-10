  
#include "ConditionalClades.hpp"
#include <trees/PLLUnrootedTree.hpp>
#include <IO/Logger.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <trees/PLLRootedTree.hpp>
#include <sstream>

struct TreeWraper {
  std::shared_ptr<PLLUnrootedTree> tree;
  corax_unode_t *root;
  size_t hash;
  double ll;
  TreeWraper(const std::string newickStr,
      bool rooted):root(nullptr), ll(0.0) {
    tree = std::make_shared<PLLUnrootedTree>(newickStr, false);
    if (rooted) {
      PLLRootedTree rootedTree(newickStr, false);
      root = tree->getVirtualRoot(rootedTree);
      hash = tree->getRootedTreeHash(root);
    } else {
      hash = tree->getUnrootedTreeHash();
    }
  }

  bool operator ==(const TreeWraper &other) const
  {
    return root == other.root && 
      PLLUnrootedTree::areIsomorphic(*tree, *(other.tree));
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
          return tree.hash;
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
    double count,
    bool useLikelihoods)
{
  assert(count);
  auto cladeCountsIt = cladeCounts.find(cid);
  if (cladeCountsIt == cladeCounts.end()) {
    cladeCounts.insert({cid, count});
  } else {
    if (useLikelihoods) {
      cladeCountsIt->second = std::max(count, cladeCountsIt->second);
    } else {
      cladeCountsIt->second += count;
    }
  }
}

static void addSubclade(CID cid,
    CID subcladeCID,
    SubcladeCounts &subcladeCounts,
    double count,
    bool useLikelihoods)
{
  addClade(subcladeCID, subcladeCounts[cid], count, useLikelihoods);
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


static void readTrees(const std::string &inputFile,
    const std::string &likelihoodFile,
    bool rooted,
    WeightedTrees &weightedTrees,
    unsigned int &inputTrees,
    unsigned int &uniqueInputTrees)
{
  std::ifstream infile(inputFile);
  std::string line;
  bool useLikelihoods = likelihoodFile.size() != 0;
  std::vector<std::string> lines;
  std::vector<double> likelihoods;
  while (std::getline(infile, line)) {
    lines.push_back(line);
  }
  if (useLikelihoods) {
    std::ifstream llFile(likelihoodFile);
    while (std::getline(llFile, line)) {
      double ll = std::stod(line);
      likelihoods.push_back(ll);
    }
    assert(likelihoods.size() == lines.size());
  }
  
  unsigned int llIndex = 0;
  for (auto line: lines) {
    TreeWraper wraper(line, rooted);
    if (useLikelihoods) {
      wraper.ll = likelihoods[llIndex++];
    };
    auto it = weightedTrees.find(wraper);
    if (it == weightedTrees.end()) {
      weightedTrees.insert({wraper, 1});
    } else {
      if (useLikelihoods) {
        if (wraper.ll > it->first.ll) {
          weightedTrees.erase(it);
          weightedTrees.insert({wraper, 1});
        }
      } else{
        it->second++;
      }
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
  CCPClade emptyClade(anyTree.getLeafNumber(), false);
  CCPClade fullClade(anyTree.getLeafNumber(), true);
  auto leafNumber = anyTree.getLeafNumber();
  std::vector<corax_unode_t *> nodes;
  std::vector<char> buffer;
  OrderedClades orderedClades;
  // add the root clade (containing all taxa)
  orderedClades.insert(fullClade);
  // iterate over all trees of the tree distribution
  for (auto pair: weightedTrees) {
    auto virtualRoot = pair.first.root;
    auto &tree = *(pair.first.tree);
    std::vector<CCPClade> nodeIndexToClade(tree.getDirectedNodeNumber(), 
        emptyClade);
    if (virtualRoot) {
      postOrderNodes.emplace_back(tree.getPostOrderNodesRooted(virtualRoot));
    } else {
      postOrderNodes.emplace_back(tree.getPostOrderNodes());
    }
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
// if CIDToDeviation is not null, fill it with the map deviations
static void fillCladeCounts(const WeightedTrees &weightedTrees,
    const CladeToCID &cladeToCID,
    const std::unordered_map<std::string, unsigned int> &leafToId,
    const std::vector<std::vector<corax_unode_t*> > &postOrderNodes,
    SubcladeCounts &subcladeCounts,
    bool useLikelihoods,
    std::unordered_map<unsigned int, double> *CIDToDeviation = nullptr)
{
  auto &anyTree = *(weightedTrees.begin()->first.tree);
  auto emptyClade = CCPClade(anyTree.getLeafNumber(), false);
  auto fullClade = CCPClade(anyTree.getLeafNumber(), true);
  auto fullCladeCID = cladeToCID.at(fullClade);
  std::vector<CCPClade> nodeIndexToClade(anyTree.getDirectedNodeNumber(), emptyClade); 
  std::vector<CID> nodeIndexToCID(anyTree.getDirectedNodeNumber());
  // root clade
  unsigned int weightedTreeIndex = 0;
  for (auto pair: weightedTrees) {
    std::vector<double> deviations;
    if (CIDToDeviation) {
      deviations = pair.first.tree->getMADRelativeDeviations();
    }
    double frequency = useLikelihoods ? pair.first.ll : pair.second;
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
          addSubclade(cid, leftCID, subcladeCounts, frequency, useLikelihoods);
        } else {
          auto rightCID = cladeToCID.at(rightClade);
          addSubclade(cid, rightCID, subcladeCounts, frequency, useLikelihoods);
        }
      }

      nodeIndexToCID[node->node_index] = cid;
      if (pair.first.root == nullptr || node == pair.first.root || node->back == pair.first.root) {
        addSubclade(fullCladeCID, cid, subcladeCounts, frequency, useLikelihoods);
      }
      if (CIDToDeviation) {
        CIDToDeviation->insert({cid, deviations[node->node_index]});
      }
    }
    weightedTreeIndex++;
  }

}




ConditionalClades::ConditionalClades(const std::string &inputFile,
    const std::string &likelihoods,
      CCPRooting ccpRooting):
  _inputTrees(0),
  _uniqueInputTrees(0),
  _ccpRooting(ccpRooting)
{
  std::ifstream is(inputFile);
  if (is.peek() == '#') {
    try {
      buildFromALEFormat(inputFile, ccpRooting);
      //printContent();
      return;
    } catch (...) {
      // maybe that wasn't the ALE format after all!
      // try reading the file as a list of gene trees
    }
  }
  buildFromGeneTrees(inputFile, likelihoods, ccpRooting);
}

static CID getCIDFromALE(size_t readCID)
{
  return readCID - 1; 
}
 
void ConditionalClades::buildFromALEFormat(const std::string &inputFile, 
        CCPRooting ccpRooting)
{
  assert(ccpRooting == CCPRooting::UNIFORM);
  std::ifstream is(inputFile);
  std::string line;
  // read comment "#constructor string"
  std::getline(is, line);
  std::string constructorString;
  std::getline(is, constructorString);
  PLLUnrootedTree tree(constructorString, false);
  auto N = tree.getLeafNumber();
  // read comment "#observations"
  std::getline(is, line);
  is >> _inputTrees; 
  std::getline(is, line);
  // read comment "#Bip counts"
  std::getline(is, line);
  // read the bip counts
  std::unordered_map<size_t, double> bipCounts; 
  size_t cladeNumber = 0;
  while (is.peek() != '#') {
    size_t id;
    size_t count;
    is >> id >> count;
    std::getline(is, line);
    cladeNumber = std::max(cladeNumber, id);
    id = getCIDFromALE(id);
    bipCounts.insert({id, count});
  }
  for (unsigned int i = 0; i < N; ++i) {
    bipCounts.insert({i, _inputTrees});
  }
  cladeNumber++; // for the clade with all taxa
  // read comment "#Bip _bls"
  std::getline(is, line);
  // skip the bls for now 
  while (is.peek() != '#') {
    std::getline(is, line);
  }
  // read comment "#Dip _counts"
  std::getline(is, line);
  // read the dips
  _allCladeSplits.resize(cladeNumber);
  while (is.peek() != '#') {
    CladeSplit split;
    std::getline(is, line);
    std::istringstream iss(line);
    iss >> split.parent >> split.left >> split.right >> split.frequency;
    split.parent = getCIDFromALE(split.parent);
    split.left = getCIDFromALE(split.left);
    split.right = getCIDFromALE(split.right);
    split.frequency /= bipCounts[split.parent];
    _allCladeSplits[split.parent].push_back(split);
  }
  // read comment "#last_leafset_id"
  std::getline(is, line);
  size_t cladeNumberCheck;
  is >> cladeNumberCheck;
  cladeNumberCheck++;
  assert(cladeNumber == cladeNumberCheck);
  std::getline(is, line);
  
  // read comment "#leaf-id"
  std::getline(is, line);
  _idToLeaf.resize(N);
  while (is.peek() != '#') {
    std::getline(is, line);
    std::istringstream iss(line);
    std::string leaf;
    size_t id;
    iss >> leaf >> id;
    id--;
    _idToLeaf[id] = leaf;
  }
  // read comment "#set-id"
  std::getline(is, line);
  _CIDToClade.resize(cladeNumber);
  std::vector<CID> cidMapping(cladeNumber);
  cidMapping[cladeNumber - 1] = cladeNumber - 1;
  CID newCID = 0;
  while (is.peek() != '#') {
    std::getline(is, line);
    std::istringstream iss(line);
    CCPClade clade(N);
    size_t CID;
    size_t leafID;
    std::string separator;
    iss >> CID;
    CID = getCIDFromALE(CID);
    cidMapping[CID] = newCID++;
    iss >> separator;
    unsigned int cladeSize = 0;
    while (!iss.eof()) {
      iss >> leafID;
      leafID--;
      clade.set(leafID);
      cladeSize++;
    }
    _CIDToClade[CID] = clade;
    _cladeToCID.insert({clade, CID});
    if (cladeSize == 1) {
      // leaf
      _CIDToLeaf[CID] = _idToLeaf[leafID];
    }
  }
  // now fill the root clade
  CCPClade fullClade(N, true);
  size_t fullCID = cladeNumber -1;
  _CIDToClade[fullCID] = fullClade;
  _cladeToCID.insert({fullClade, fullCID});
  for (unsigned int cid = 0; cid < _CIDToClade.size() - 1; ++cid) {
    CladeSplit split;
    split.parent = fullCID;
    auto clade = _CIDToClade[cid];
    if (!clade[0]) {
      // do not add both a split and its complementary
    //  continue;
    }
    auto compClade = ~clade;
    split.left = cid;
    split.right = _cladeToCID[compClade];
    assert(bipCounts[cid] == bipCounts[split.right]);
    split.frequency = double(bipCounts[cid]) / double(_inputTrees * 2 * (2 * N - 3));
    _allCladeSplits[split.parent].push_back(split);
  }
  for (auto &splits: _allCladeSplits) {
    if (!splits.size()) {
      continue;
    }
    double sum = 0.0;
    for (auto &split: splits) {
      sum += split.frequency;
    }
    assert(fabs(sum - 1.0) < 0.0000001);
  }
  reorderClades(cidMapping);
}

void ConditionalClades::buildFromGeneTrees(const std::string &inputFile,
    const std::string &likelihoods,
    CCPRooting ccpRooting)
{
  std::unordered_map<std::string, unsigned int> leafToId;
  CCPClade emptyClade;
  CCPClade fullClade;
  WeightedTrees weightedTrees;
  readTrees(inputFile,
      likelihoods,
      (ccpRooting == CCPRooting::ROOTED),
      weightedTrees, 
      _inputTrees, 
      _uniqueInputTrees); 
  auto &anyTree = *(weightedTrees.begin()->first.tree);
  for (auto leaf: anyTree.getLabels()) {
    leafToId.insert({leaf, leafToId.size()});
    _idToLeaf.push_back(leaf);
  }
  emptyClade = CCPClade(anyTree.getLeafNumber(), false);
  fullClade = CCPClade(anyTree.getLeafNumber(), true);
  
  std::vector<std::vector<corax_unode_t*> > postOrderNodes;
  // first pass to get all clades and assign them a CID
  extractClades(weightedTrees, 
      leafToId,
      _cladeToCID,
      _CIDToClade,
      _CIDToLeaf,
      postOrderNodes
      );


  SubcladeCounts subcladeCounts(_CIDToClade.size());
  std::unique_ptr<std::unordered_map<unsigned int, double> >
    CIDToDeviation;
  if (madRooting()) {
    CIDToDeviation = 
      std::make_unique<std::unordered_map<unsigned int, double> >();
    }
  // second pass to count the number of occurence of each clade,
  // and for each clade, the number of occurences of its subclades
  bool useLikelihoods = likelihoods.size();
  fillCladeCounts(weightedTrees,
      _cladeToCID,
      leafToId,
      postOrderNodes,
      subcladeCounts,
      useLikelihoods,
      CIDToDeviation.get());
  _fillCCP(subcladeCounts, useLikelihoods, CIDToDeviation.get());
}
 

void normalizeFrequencies(CladeSplits &splits, bool logScale) {
  
  if (logScale) {
    // rescale to loglikelihoods to likelihoods
    double max = -std::numeric_limits<double>::max();
    for (const auto &split: splits) {
      max = std::max(max, split.frequency);
    }
    for (auto &split: splits) {
      split.frequency = std::max(0.000001, exp(split.frequency - max));
      //Logger::info << split.frequency << " ";
    }
    //Logger::info << std::endl;
  } 
  double sum = 0.0;
  for (const auto &split: splits) {
    sum += split.frequency;
  }
  for (auto &split: splits) {
    split.frequency /= sum;
  }
}
  
void ConditionalClades::_fillCCP(SubcladeCounts &subcladeCounts,
    bool useLikelihoods,
    std::unordered_map<unsigned int, double> *CIDToDeviation)
{
  _allCladeSplits.clear();
  auto cladesNumber = _CIDToClade.size();
  _allCladeSplits.resize(cladesNumber);
  auto rootCID = cladesNumber -  1;
  CCPClade fullClade(_CIDToClade[0].size(), true);
  assert(rootCID = _cladeToCID[fullClade]);
  for (unsigned int cid = 0; cid < cladesNumber; ++cid) {
    auto &clade = _CIDToClade[cid];
    auto &cladeSplits = _allCladeSplits[cid];
    if (subcladeCounts[cid].size()) {
      // internal clade
      for (auto &subcladeCount: subcladeCounts[cid]) {
        auto CIDLeft = subcladeCount.first;
        auto cladeLeft = _CIDToClade[CIDLeft];
        auto cladeRight = getComplementary(clade, cladeLeft);
        if (cladeRight < cladeLeft) {
          continue;
        }
        auto CIDRight = _cladeToCID[cladeRight];
        CladeSplit split;
        split.parent = cid;
        split.left = CIDLeft;
        split.right = CIDRight;
        double frequency = double(subcladeCount.second); 
        if (CIDToDeviation && split.parent == rootCID) {
          double deviation = (*CIDToDeviation)[split.left] + 1.0;
          frequency /= (deviation * deviation);
        }
        split.frequency = frequency;
        cladeSplits.push_back(split);
      }
      normalizeFrequencies(cladeSplits, useLikelihoods);
    }
  }
}

void ConditionalClades::printContent() const
{

  for (unsigned int CID = 0; CID < _allCladeSplits.size(); ++CID) {
    auto &clade = _CIDToClade[CID];
    auto &cladeSplits = _allCladeSplits[CID];
    std::cerr << "clade " << CID << " ";
    printClade(clade, _idToLeaf);
    std::cerr << std::endl;
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
  //printContent();
}

void ConditionalClades::printStats() const
{
  std::cerr << "Leaves " << _CIDToLeaf.size() << " clades: " << _cladeToCID.size() <<  " unique trees " << _uniqueInputTrees << std::endl;

}

void ConditionalClades::reorderClades(const std::vector<CID> &mappings)
{
  auto cladeNumber = mappings.size();
  // remap _CIDToLeaf
  CIDToLeaf newCIDToLeaf(cladeNumber);
  for (const auto &pair: _CIDToLeaf) {
    newCIDToLeaf[mappings[pair.first]] = pair.second;
  }
  _CIDToLeaf = newCIDToLeaf;
  // remap _CIDToClade
  for (const auto &pair: _cladeToCID) {
    _CIDToClade[mappings[pair.second]] = pair.first;
  }
  // remap _cladeToCID
  for (CID cid = 0; cid < cladeNumber; ++cid) {
    _cladeToCID.insert({_CIDToClade[cid], cid});
  }
  // remap _allCladeSplits
  std::vector<CladeSplits> newAllCladeSplits(cladeNumber);
  for (unsigned int oldCid = 0; oldCid < cladeNumber; ++oldCid) {
    auto newCid = mappings[oldCid];
    const auto &cladeSplits = _allCladeSplits[oldCid];
    for (const auto &split: cladeSplits) {
      CladeSplit newSplit;
      newSplit.parent = mappings[split.parent];
      newSplit.left = mappings[split.left];
      newSplit.right = mappings[split.right];
      newSplit.frequency = split.frequency;
      newAllCladeSplits[newCid].push_back(newSplit);
    }
  }
  _allCladeSplits = newAllCladeSplits;
}
