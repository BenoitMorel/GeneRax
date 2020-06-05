#include <trees/PLLRootedTree.hpp>
#include <cassert>
#include <iostream>
#include <fstream>

const std::string tree1("((A,B),(C,D));");
const std::string tree2("(E, ((A,B),(C,D)));");


const std::string tree3("(((A,B),(C,D)),((E,(F,G)),(((H,I),J),((K,L),(((M,N),(O,P)),(Q,(R,(S,(T,U)))))))));");


void testTree(const std::string &treeStr)
{
  std::ofstream os("temp.txt");
  os << treeStr << std::endl;
  os.close();
  PLLRootedTree treeFromFile("temp.txt", true);
  PLLRootedTree tree(treeStr, false);
  assert(tree.getNodesNumber() == treeFromFile.getNodesNumber());
}




int main(int, char**)
{
  testTree(tree1); 
  testTree(tree2); 
  testTree(tree3); 
  //testTree(frediTreeStr); 
  std::cout << "Test PLLRootedTree ok!" << std::endl;
  return 0;
}






