#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>
#include <chrono>


void printLibpllNode(pll_unode_s *node, ostream &os, bool isRoot)
{
    if (node->next) {
        os << "(";
        printLibpllNode(node->next->back, os, false);
        os << ",";
        printLibpllNode(node->next->next->back, os, false);
        os << ")";
    } else {
        os << node->label;
    }
    os << ":" << (isRoot ? node->length / 2.0 : node->length);
}

void printLibpllTreeRooted(pll_unode_t *root, ostream &os){
    os << "(";
    printLibpllNode(root, os, true);
    os << ",";
    printLibpllNode(root->back, os, true);
    os << ");" << endl;
}

BPPNode createNode(BPPTree tree, BPPNode father = 0, const std::string&name = "", double branchLength = 0.0) {
  BPPNode node(new bpp::PhyloNode(name));
  if (!father) {
      tree->createNode(node);
  } else {
      BPPBranch branch(new bpp::PhyloBranch);
      if (tree->hasFather(father)) {
          branch->setLength(branchLength);
      } else {
          // special case: this branch corresponds to half the branch
          // of the virtual root in the unrooted tree
          branch->setLength(branchLength / 2.0);
      }
      tree->createNode(father, node, branch);
  }
  return node;
}

void addFromLibpll(BPPTree tree, BPPNode bppFather, pll_unode_s *libpllNode)
{
    BPPNode newNode = createNode(tree, bppFather, libpllNode->label ? libpllNode->label : "", libpllNode->length);
    if (libpllNode->next) {
        addFromLibpll(tree, newNode, libpllNode->next->back);
        addFromLibpll(tree, newNode, libpllNode->next->next->back);
    }
}

std::shared_ptr<bpp::PhyloTree> buildFromLibpll(std::shared_ptr<LibpllEvaluation> evaluation, pll_unode_s *libpllRoot)
{
    std::shared_ptr<bpp::PhyloTree> tree(new bpp::PhyloTree(true));
    BPPNode root = createNode(tree, 0, "root");
    tree->setRoot(root);
    addFromLibpll(tree, root, libpllRoot);
    addFromLibpll(tree, root, libpllRoot->back);
    return tree;
}
