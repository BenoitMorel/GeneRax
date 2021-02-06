  
#include "ConditionalClades.hpp"
#include <trees/PLLUnrootedTree.hpp>

#include <fstream>
#include <iostream>


static CCPClade getComplementary(const CCPClade &clade,
    const CCPClade &subclade)
{
  CCPClade complementary(clade);
  for (unsigned int i = 0; i < clade.size(); ++i) {
    complementary[i] = clade[i] && !subclade[i];
  }
  return complementary;
}

static void addClade(const CCPClade &clade,
    CladeCounts &cladeCounts)
{
  auto cladeCountsIt = cladeCounts.find(clade);
  if (cladeCountsIt == cladeCounts.end()) {
    cladeCounts.insert({clade, 1});
  } else {
    cladeCountsIt->second++;
  }
}

static void addSubclade(const CCPClade &clade,
    const CCPClade &subclade,
    SubcladeCounts &subcladeCounts)
{
  auto subcladeCountsIt = subcladeCounts.find(clade);
  if (subcladeCountsIt == subcladeCounts.end()) {
    CladeCounts cladeCounts;
    subcladeCounts.insert({clade, CladeCounts()});
    subcladeCountsIt = subcladeCounts.find(clade);
  }
  addClade(subclade, subcladeCountsIt->second);
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
  while (std::getline(infile, line)) {
    // todo: support empty lines
    // tododetect duplicated trees
    PLLUnrootedTree tree(line, false);
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
        clade[id] = true; 
      } else {
        auto &leftClade = nodeIndexToClade[node->next->back->node_index];
        auto &rightClade = nodeIndexToClade[node->next->next->back->node_index];
        // todo: WE NEED A PROPER BITSET!!!
        for (unsigned int i = 0; i < clade.size(); ++i) {
          clade[i] = leftClade[i] || rightClade[i];
        }
        if (leftClade > rightClade) {
          addSubclade(clade, leftClade, subcladeCounts);
        } else {
          addSubclade(clade, rightClade, subcladeCounts);
        }
      }
      orderedClades.insert(clade);

      addClade(clade, cladeCounts);
      
      // root clade
      addClade(fullClade, cladeCounts);
      addSubclade(fullClade, clade, subcladeCounts);
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

