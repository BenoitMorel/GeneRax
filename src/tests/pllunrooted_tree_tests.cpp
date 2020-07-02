#include <trees/PLLUnrootedTree.hpp>
#include <cassert>
#include <iostream>
#include <fstream>
#include <unordered_map>

const std::string tree1("((A,B),(C,D), (E,F));");
const std::string tree2("(E, F, ((A,B),(C,D)));");

const std::string treeWithDistances1 =
  "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.2,(E:0.2,F:0.3):0.1):0.0;";

void testTree(const std::string &treeStr)
{
  std::ofstream os("temp.txt");
  os << treeStr << std::endl;
  os.close();
  PLLUnrootedTree treeFromFile("temp.txt", true);
  PLLUnrootedTree tree(treeStr, false);
  assert(tree.getNodesNumber() == treeFromFile.getNodesNumber());
}


bool almostEqual(double x, double y) 
{
  return fabs(x - y) < 0.000000001;
}

void testDistances()
{
  PLLUnrootedTree t1(treeWithDistances1, false);
  MatrixDouble distances;
  t1.computePairwiseDistances(distances);
  std::unordered_map<std::string, int> indices;
  for (auto leaf: t1.getLeaves()) {
    indices.insert({std::string(leaf->label), leaf->node_index});
  }
  assert(indices.size() == t1.getLeavesNumber());
  
  for (unsigned int i = 0; i < t1.getLeavesNumber(); ++i) {
    assert(distances[i][i] == 0.0);
    for (unsigned int j = 0; j < t1.getLeavesNumber(); ++j) {
      assert(almostEqual(distances[i][j], distances[j][i]));
    }
  }
  assert(almostEqual(distances[indices["A"]][indices["B"]],  0.2));
  assert(almostEqual(distances[indices["A"]][indices["C"]],  0.7));
  assert(almostEqual(distances[indices["A"]][indices["E"]],  0.6));
  assert(almostEqual(distances[indices["A"]][indices["F"]],  0.7));
}



int main(int, char**)
{
  testTree(tree1); 
  testTree(tree2); 
  testTree(treeWithDistances1);
  testDistances();
  std::cout << "Test PLLUnrootedTree ok!" << std::endl;
  return 0;
}







