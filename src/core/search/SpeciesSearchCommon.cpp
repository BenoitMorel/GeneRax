#include "SpeciesSearchCommon.hpp"

#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>

static std::string getSubtreeID(pll_rnode_t *subtree)
{
  if (!subtree->left) {
    return std::string(subtree->label);
  }
  std::string res("(");
  std::string id1 = getSubtreeID(subtree->left);
  std::string id2 = getSubtreeID(subtree->right);
  if  (id1 > id2) {
    std::swap(id1, id2);
  }
  return std::string("(") + id1 + "," + id2 + ")";
}

void RootLikelihoods::saveValue(pll_rnode_t *subtree, double ll) 
{
  auto id = getSubtreeID(subtree); 
  idToLL[id] = ll;
}

void RootLikelihoods::fillTree(PLLRootedTree &tree)
{
  std::vector<double> nodeIdToLL(tree.getNodesNumber(), 0.0);
  double bestLL = -std::numeric_limits<double>::infinity();
  for (auto node: tree.getNodes()) {
    auto id = getSubtreeID(node);
    if (idToLL.find(id) != idToLL.end()) {
      // we have a likelihood value
      auto value = idToLL[id];
      nodeIdToLL[node->node_index] = value;
      bestLL = std::max<double>(value, bestLL);
    }
  }
  for (auto node: tree.getNodes()) {
    bool hasValue = nodeIdToLL[node->node_index] != 0.0;
    std::string label;
    if (hasValue) {
      double value = nodeIdToLL[node->node_index] - bestLL;
      if(!node->left && node->label) {
        label = std::string(node->label);
        free(node->label);
        node->label = nullptr;
        label += "_";
      }
      label += std::to_string(value);
      node->label = (char*)malloc(label.size() + 1);
      memcpy(node->label, label.c_str(), label.size());
      node->label[label.size()] = 0;
    }
  }
}

