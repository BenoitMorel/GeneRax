  
#include "ConditionalClades.hpp"
#include <trees/PLLUnrootedTree.hpp>

#include <fstream>
#include <iostream>




static CCPClade getComplementary(const CCPClade &clade,
    const CCPClade &subclade)
{
#ifdef GENESISIMPLEM
  CCPClade complementary(subclade);
  complementary.negate();
  complementary &= clade;
#else
  CCPClade complementary(clade);
  for (unsigned int i = 0; i < clade.size(); ++i) {
    complementary[i] = clade[i] && !subclade[i];
  }
#endif
  return complementary;
}

static void addClade(const CCPClade &clade,
    CladeCounts &cladeCounts,
    unsigned int count)
{
  auto cladeCountsIt = cladeCounts.find(clade);
  if (cladeCountsIt == cladeCounts.end()) {
    cladeCounts.insert({clade, count});
  } else {
    cladeCountsIt->second += count;
  }
}

static void addSubclade(const CCPClade &clade,
    const CCPClade &subclade,
    SubcladeCounts &subcladeCounts,
    unsigned int count)
{
  auto subcladeCountsIt = subcladeCounts.find(clade);
  if (subcladeCountsIt == subcladeCounts.end()) {
    CladeCounts cladeCounts;
    subcladeCounts.insert({clade, CladeCounts()});
    subcladeCountsIt = subcladeCounts.find(clade);
  }
  addClade(subclade, subcladeCountsIt->second, count);
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

struct WeightedTree {
  WeightedTree(std::shared_ptr<PLLUnrootedTree> tree, 
      unsigned int count): tree(tree), count(count) {}

  std::shared_ptr<PLLUnrootedTree> tree;
  unsigned int count;
};


ConditionalClades::ConditionalClades(const std::string &newickFile)
{
  std::ifstream infile(newickFile);
  std::string line;
  std::unordered_map<std::string, unsigned int> leafToId;
  CCPClade emptyClade;
  CCPClade fullClade;
  CladeCounts cladeCounts;
  SubcladeCounts subcladeCounts;
  OrderedClades orderedClades;
  std::unordered_map<size_t, WeightedTree> weightedTrees;
  while (std::getline(infile, line)) {
    // todo: support empty lines
#undef CCPCOMPRESS 
    // BEFORE ENABLING CCPCOMPRESS: check no collision
    // with hash or add an == operator for unrooted trees
#ifdef CCPCOMPRESS
    auto tree = std::make_shared<PLLUnrootedTree>(line, false);
    auto hash = tree->getUnrootedTreeHash();
    if (weightedTrees.find(hash) == weightedTrees.end()) {
      WeightedTree weightedTree(tree, 1);
      weightedTrees.insert({hash, weightedTree});
    } else {
      weightedTrees.at(hash).count++;
    }
  }
  std::cout << "Number of different trees: " << weightedTrees.size() << std::endl;
  for (auto pair: weightedTrees) {
    auto &tree = *(pair.second.tree);
    auto treeCount = pair.second.count;
#else
    PLLUnrootedTree tree(line, false);
    unsigned int treeCount = 1;
#endif
    if (leafToId.size() == 0) {
      // first tree, initiate mappings
      for (auto leaf: tree.getLabels()) {
        leafToId.insert({leaf, leafToId.size()});
        _idToLeaf.push_back(leaf);
      }
      emptyClade = CCPClade(tree.getLeavesNumber(), false);
      fullClade = CCPClade(tree.getLeavesNumber(), true);
      orderedClades.insert(fullClade);
    }
    std::vector<CCPClade> nodeIndexToClade(tree.getDirectedNodesNumber(), 
        emptyClade);
    for (auto node: tree.getPostOrderNodes()) {
      auto nodeIndex = node->node_index;
      auto &clade = nodeIndexToClade[nodeIndex];
      if (!node->next) { // leaf
        auto id = leafToId[node->label];
#ifdef GENESISIMPLEM
        clade.set(id);
#else
        clade[id] = true; 
#endif
      } else {
        auto &leftClade = nodeIndexToClade[node->next->back->node_index];
        auto &rightClade = nodeIndexToClade[node->next->next->back->node_index];
#ifdef GENESISIMPLEM
        clade = leftClade | rightClade;
#else
        for (unsigned int i = 0; i < clade.size(); ++i) {
          clade[i] = leftClade[i] || rightClade[i];
        }
#endif
        if (leftClade < rightClade) {
          addSubclade(clade, leftClade, subcladeCounts, treeCount);
        } else {
          addSubclade(clade, rightClade, subcladeCounts, treeCount);
        }
      }
      orderedClades.insert(clade);

      addClade(clade, cladeCounts, treeCount);
      
      // root clade
      addClade(fullClade, cladeCounts, treeCount);
      // todo should we check that a clade/complementary is only
      // added once to the root??
      addSubclade(fullClade, clade, subcladeCounts, treeCount);
    }
  }
  _fillCCP(cladeCounts, subcladeCounts, orderedClades);
}
  
void ConditionalClades::_fillCCP(CladeCounts &cladeCounts,
      SubcladeCounts &subcladeCounts,
      OrderedClades &orderedClades)
{
  _allCladeSplits.clear();
  _allCladeSplits.resize(orderedClades.size());
  for (auto it = orderedClades.begin(); it != orderedClades.end(); ++it) {
    auto &clade = *it;
    auto cladeCount = cladeCounts[clade];
    unsigned int CID = _CIDToClade.size();
    auto &cladeSplits = _allCladeSplits[CID];
    _CIDToClade.push_back(clade);
    _cladeToCID[clade] = CID;
    auto subcladeCountIt = subcladeCounts.find(clade);
    if (subcladeCountIt != subcladeCounts.end()) {
      // internal clade
      double sumFrequencies = 0.0;
      for (auto &subcladeCount: subcladeCountIt->second) {
        auto &cladeLeft = subcladeCount.first;
        auto cladeRight = getComplementary(clade, cladeLeft);
        auto CIDLeft = _cladeToCID[cladeLeft];
        auto CIDRight = _cladeToCID[cladeRight];
        double frequency = double(subcladeCount.second) / double(cladeCount); 
        sumFrequencies += frequency;
        CladeSplit split;
        split.parent = CID;
        split.left = CIDLeft;
        split.right = CIDRight;
        split.frequency = frequency;
        cladeSplits.push_back(split);
      }
      // Check that frequencies sum to one
      assert(fabs(1.0 - sumFrequencies) < 0.000001);
    } else {
      // leaf clade
      unsigned int pos = 0;
      // todo: make it faster with a bitset
      for (pos = 0; !clade[pos]; ++pos) {}
      _CIDToLeaf[CID] = _idToLeaf[pos];
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

