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
  
void TreeDuplicatesFinder::findDuplicates(PerCoreGeneTrees &perCoreGeneTrees)
{

  std::unordered_map<pll_unode_t *, unsigned int> nodeClassIdentifiers; 
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
      if (nodeClassIdentifiers.find(node) != nodeClassIdentifiers.end()) {
        // we already processed this node
        continue;
      }
      std::string nodeClassString;
      if (!node->next) {
        nodeClassString = leafToSpecies[i];     
      } else {
        auto sonLeft = node->next->back;
        auto sonRight = node->next->next->back;
        unsigned int sonLeftId = nodeClassIdentifiers[sonLeft];
        unsigned int sonRightId = nodeClassIdentifiers[sonRight];
        // ordering does not matter for equality
        if (sonLeftId > sonRightId) {
          std::swap(sonLeftId, sonRightId);
        }
        nodeClassString = std::to_string(sonLeftId) + "-" + std::to_string(sonRightId);
      }
      if (stringToIdentifier.find(nodeClassString) == stringToIdentifier.end()) {
        stringToIdentifier[nodeClassString] = newClassIdentifier++;
      }
      nodeClassIdentifiers[node] = stringToIdentifier[nodeClassString];
      if (hasParent(node)) {
        nextDepthNodes.push_back(getParentLeft(node));
        nextDepthNodes.push_back(getParentRight(node));
      }
    }
    // compute next depth nodes
    std::swap(currentDepthNodes, nextDepthNodes);
  }
  Logger::info << "Subtrees: " << nodeClassIdentifiers.size() << std::endl;
  Logger::info << "Unique subtrees: " << newClassIdentifier << std::endl;
}


