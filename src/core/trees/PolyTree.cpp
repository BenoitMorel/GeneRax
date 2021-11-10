#include "PolyTree.hpp"


#include <iostream>

static bool contractBranch(corax_unode_t *branch)
{
  return branch->length <= 0.000001 && branch->next;
}

static void fillChildrenRec(corax_unode_t *node,
    std::vector<corax_unode_t*> &children)
{
  auto left = node->next->back;
  auto right = node->next->next->back;
  assert(left);
  assert(right);
  if (!contractBranch(left)) {
    children.push_back(left);
  } else {
    fillChildrenRec(left, children);
  }
  if (!contractBranch(right)) {
    children.push_back(right);
  } else {
    fillChildrenRec(right, children);
  }

}

PolyTree::PolyTree(PLLUnrootedTree &tree):
  _tree(tree),
  _cells(tree.getDirectedNodesNumber())
{
  for (auto node: tree.getPostOrderNodes()) {
    auto gid = node->node_index;
    auto &cell = _cells[gid];
    cell.gid = gid;
    if (node->next) {
      fillChildrenRec(node, cell.children);
    }
  }
}


static std::string buildStringRec(corax_unode_t *node,
    const std::vector<PolyTree::Cell> &cells)
{
  auto gid = node->node_index;
  std::string res;
  auto &cell = cells[gid];
  if (cell.children.size()) {
    res += "(";
    for (unsigned int i = 0; i < cell.children.size(); ++i) {
      res += buildStringRec(cell.children[i], cells);
      if (i != cell.children.size() - 1) {
        res += ",";
      }
    }
    res += ")";
  }
  if (node->label) {
    res += std::string(node->label);
  }
  res += ":";
  res += std::to_string(node->length);
  return res;
}

std::string PolyTree::getNewickString() const
{
  auto root = _tree.getNode(0);
  auto res1 = buildStringRec(root, 
      _cells);
  auto res2 = buildStringRec(root->back,
      _cells);
  auto res = "(" + res1 + "," + res2 + ");";
  return res;
}

