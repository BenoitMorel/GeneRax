  
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
    std::vector<std::string> &idToLeaf)
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
    }
    std::cerr << "Treating tree " << tree.getNewickString() << std::endl;
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
          addSubclade(clade, leftClade, _subcladeCounts);
        } else {
          addSubclade(clade, rightClade, _subcladeCounts);
        }
        _internalClades.insert(clade);
      }
      addClade(clade, _cladeCounts);
    }
  }
  printContent();
}

void ConditionalClades::printContent()
{
  for (auto it = begin(); it != end(); ++it) {
    auto &clade = *it;
    auto &localCladeCounts = _subcladeCounts[clade];
    printClade(clade, _idToLeaf);
    std::cerr << std::endl;
    std::cerr << "Count = " << _cladeCounts[clade] << std::endl;
    for (auto &localCladeCount: localCladeCounts) {
      std::cerr << "  ";
      auto &cladeLeft = localCladeCount.first;
      auto cladeRight = getComplementary(clade, cladeLeft);
      printClade(cladeLeft, _idToLeaf);
      std::cerr << "|";
      printClade(cladeRight, _idToLeaf);
      std::cerr << " -> " << localCladeCount.second << std::endl;
    }
    std::cerr << std::endl;
  }
}

