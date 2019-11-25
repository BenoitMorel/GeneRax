#include "TreeDuplicatesFinder.hpp"
#include <trees/PerCoreGeneTrees.hpp>
#include <IO/Logger.hpp>

static bool hasParent(pll_unode_t *node) 
{
  return node->back->next;
}

static pll_unode_t *getParentLeft(pll_unode_t *node)
{
  return node->back->next;
}

static pll_unode_t *getParentRight(pll_unode_t *node) 
{
  return node->back->next->next;
}
  
void TreeDuplicatesFinder::findDuplicates(PerCoreGeneTrees &perCoreGeneTrees, TreeDuplicates &duplicates)
{
  std::vector<pll_unode_t *> currentDepthNodes;
  std::vector <std::string> leafToSpecies;
  unsigned int newClassIdentifier = 1;
  Logger::info << "Families: " << perCoreGeneTrees.getTrees().size() << std::endl;
  // fill leaves
  for (auto &geneTree: perCoreGeneTrees.getTrees()) {
    auto tree = geneTree.geneTree;
    for (auto leaf: tree->getLeaves()) {
      currentDepthNodes.push_back(leaf);
      leafToSpecies.push_back(geneTree.mapping.getSpecies(leaf->label));
    }
  }
  Logger::info << "Total number of leaves: " << currentDepthNodes.size() << std::endl;
  while (currentDepthNodes.size()) {
    std::vector<pll_unode_t *> nextDepthNodes;
    std::unordered_map<std::string, unsigned int> stringToIdentifier;
    // compute their identifiers
    for (unsigned int i = 0; i < currentDepthNodes.size(); ++i) {
      auto node = currentDepthNodes[i];
      if (duplicates.find(node) != duplicates.end()) {
        // we already processed this node
        continue;
      }
      std::string nodeClassString;
      if (!node->next) {
        nodeClassString = leafToSpecies[i];     
      } else {
        auto sonLeft = node->next->back;
        auto sonRight = node->next->next->back;
        unsigned int sonLeftId = duplicates[sonLeft];
        unsigned int sonRightId = duplicates[sonRight];
        // ordering does not matter for equality
        if (sonLeftId > sonRightId) {
          std::swap(sonLeftId, sonRightId);
        }
        nodeClassString = std::to_string(sonLeftId) + "-" + std::to_string(sonRightId);
      }
      if (stringToIdentifier.find(nodeClassString) == stringToIdentifier.end()) {
        stringToIdentifier[nodeClassString] = newClassIdentifier++;
      }
      duplicates[node] = stringToIdentifier[nodeClassString];
      if (hasParent(node)) {
        nextDepthNodes.push_back(getParentLeft(node));
        nextDepthNodes.push_back(getParentRight(node));
      }
    }
    std::swap(currentDepthNodes, nextDepthNodes);
  }
  Logger::info << "Subtrees: " << duplicates.size() << std::endl;
  Logger::info << "Unique subtrees: " << newClassIdentifier << std::endl;
}

template<typename Value>
SubtreeCache<Value>::SubtreeCache(const TreeDuplicates &duplicates):
  _duplicates(duplicates),
  _values(duplicates.size()),
  _isValid(duplicates.size(), false)
{

}

template<typename Value>
void SubtreeCache<Value>::resetAll()
{
  std::fill(_isValid.begin(), _isValid.end(), false);
}

