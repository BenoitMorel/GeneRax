#include <trees/PLLRootedTree.hpp>
#include <cassert>
#include <iostream>
#include <fstream>

const std::string tree1("((A,B),(C,D));");
const std::string tree2("(E, ((A,B),(C,D)));");


const std::string tree3("(((A,B),(C,D)),((E,(F,G)),(((H,I),J),((K,L),(((M,N),(O,P)),(Q,(R,(S,(T,U)))))))));");

pll_rnode_t * getLCA(pll_rnode_t *n1,
    pll_rnode_t *n2)
{
  std::unordered_set<pll_rnode_t *> n1Ancestors;
  while (n1) {
    n1Ancestors.insert(n1);
    n1 = n1->parent;
  }
  while (n2) {
    if (n1Ancestors.find(n2) != n1Ancestors.end()) {
      return n2;
    }
    n2 = n2->parent;
  }
  assert(false);
}

void testTree(const std::string &treeStr)
{
  std::ofstream os("temp.txt");
  os << treeStr << std::endl;
  os.close();
  PLLRootedTree treeFromFile("temp.txt", true);
  PLLRootedTree tree(treeStr, false);
  assert(tree.getNodesNumber() == treeFromFile.getNodesNumber());

  for (auto n1: tree.getNodes()) {
    for (auto n2: tree.getNodes()) {
      assert(tree.getLCA(n1, n2) == getLCA(n1, n2));
    }
  }
}


int main(int, char**)
{
  testTree(tree1); 
  testTree(tree2); 
  testTree(tree3); 
  std::cout << "Test PLLRootedTree ok!" << std::endl;
  return 0;
}






