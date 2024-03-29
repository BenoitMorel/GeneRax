#include <trees/PLLUnrootedTree.hpp>
#include <cassert>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>

const std::string tree1("((A,B),(C,D), (E,F));");
const std::string tree2("(E, F, ((A,B),(C,D)));");

const std::string treeWithDistances1 =
  "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.2,(E:0.2,F:0.3):0.1):0.0;";

void testGetClade(corax_unode_t *node, unsigned int leavesNumber)
{
  auto c1 = PLLUnrootedTree::getClade(node);
  auto c2 = PLLUnrootedTree::getClade(node->back);
  for (auto e1: c1) {
    for (auto e2: c2) {
      assert(e1 != e2);
    }
  }
  assert(c1.size() + c2.size() == leavesNumber);
}

void testTree(const std::string &treeStr)
{
  std::ofstream os("temp.txt");
  os << treeStr << std::endl;
  os.close();
  PLLUnrootedTree treeFromFile("temp.txt", true);
  PLLUnrootedTree tree(treeStr, false);
  for (auto node: tree.getNodes()) {
    testGetClade(node, tree.getLeafNumber());
  }
  assert(tree.getNodeNumber() == treeFromFile.getNodeNumber());
}


bool almostEqual(double x, double y) 
{
  return fabs(x - y) < 0.000000001;
}

void testDistances(bool leavesOnly)
{
  PLLUnrootedTree t1(treeWithDistances1, false);
  MatrixDouble distances;
  t1.computePairwiseDistances(distances, leavesOnly);
  std::unordered_map<std::string, int> indices;
  for (auto leaf: t1.getLeaves()) {
    indices.insert({std::string(leaf->label), leaf->node_index});
  }
  assert(indices.size() == t1.getLeafNumber());
  
  for (unsigned int i = 0; i < t1.getLeafNumber(); ++i) {
    assert(distances[i][i] == 0.0);
    for (unsigned int j = 0; j < t1.getLeafNumber(); ++j) {
      assert(almostEqual(distances[i][j], distances[j][i]));
    }
  }
  assert(almostEqual(distances[indices["A"]][indices["B"]],  0.2));
  assert(almostEqual(distances[indices["A"]][indices["C"]],  0.7));
  assert(almostEqual(distances[indices["A"]][indices["E"]],  0.6));
  assert(almostEqual(distances[indices["A"]][indices["F"]],  0.7));
  if (!leavesOnly) {
    auto A = t1.getNode(indices["A"]);
    auto AB1 = A->back->node_index;
    auto AB2 = A->back->next->node_index;
    auto AB3 = A->back->next->next->node_index;
    assert(almostEqual(distances[AB1][indices["A"]], 0.1));
    assert(almostEqual(distances[AB2][indices["A"]], 0.1));
    assert(almostEqual(distances[AB3][indices["A"]], 0.1));
    assert(almostEqual(distances[AB1][indices["F"]], 0.6));
    assert(almostEqual(distances[AB2][indices["F"]], 0.6));
    assert(almostEqual(distances[AB3][indices["F"]], 0.6));
  }
}

void testMADRooting()
{
  PLLUnrootedTree t1(treeWithDistances1, false);
  auto deviations = t1.getMADRelativeDeviations(); 
  assert(deviations.size() == t1.getDirectedNodeNumber());
  auto minIt = std::min_element(deviations.begin(), deviations.end());
  auto minIndex = std::distance(deviations.begin(), minIt);
  corax_unode_t *root = nullptr;
  for (auto node: t1.getNodes()) {
    if (node->node_index == minIndex) {
      root = node;
      break;
    }
    if (node->next) {
      if (node->next->node_index == minIndex) {
        root = node->next;
        break;
      } 
      if (node->next->next->node_index == minIndex) {
        root = node->next->next;
        break;
      }
    }
  }
  assert(root);
  auto c1 = PLLUnrootedTree::getClade(root);
  auto c2 = PLLUnrootedTree::getClade(root->back);
  std::unordered_map<std::string, int> indices;
  for (auto leaf: t1.getLeaves()) {
    indices.insert({std::string(leaf->label), leaf->node_index});
  }
  if (c1.find(indices["A"]) == c1.end()) {
    std::swap(c1, c2);
  }
  assert(c1.find(indices["A"]) != c1.end());
  assert(c1.find(indices["B"]) != c1.end());
  assert(c1.find(indices["E"]) != c1.end());
  assert(c1.find(indices["F"]) != c1.end());
  assert(c2.find(indices["C"]) != c2.end());
  assert(c2.find(indices["D"]) != c2.end());
}

void testIsomorphism()
{
  std::vector<std::shared_ptr<PLLUnrootedTree> > class1;
  std::vector<std::shared_ptr<PLLUnrootedTree> > unclassed;
  class1.push_back(std::make_shared<PLLUnrootedTree>(
        "((A,B),(C,D),(E,F));", false));
  class1.push_back(std::make_shared<PLLUnrootedTree>(
        "((F,E),(D,C),(A,B));", false));
  class1.push_back(std::make_shared<PLLUnrootedTree>(
        "(E, F, ((A,B),(C,D)));", false));
  unclassed.push_back(std::make_shared<PLLUnrootedTree>(
        "((A,B), (C, D), E);", false));
  unclassed.push_back(std::make_shared<PLLUnrootedTree>(
        "((A,B), (C, E), (D,F));", false));
  unclassed.push_back(std::make_shared<PLLUnrootedTree>(
        "((A,C), (B, D), (E,F));", false));
  unclassed.push_back(std::make_shared<PLLUnrootedTree>(
        "((A,B), (C, D), (E,GO));", false));
  for (auto t1: class1) {
    for (auto t2: class1) {
      if (!PLLUnrootedTree::areIsomorphic(*t1, *t2)) {
        std::cerr << "Error: trees:" << std::endl;
        std::cerr << t1->getNewickString() << std::endl;
        std::cerr << t2->getNewickString() << std::endl;
        std::cerr << "Should be isomorphic" << std::endl;
        assert(false);
      }
    }
    for (auto t2: unclassed) {
      if (PLLUnrootedTree::areIsomorphic(*t1, *t2)) {
        std::cerr << "Error: trees:" << std::endl;
        std::cerr << t1->getNewickString() << std::endl;
        std::cerr << t2->getNewickString() << std::endl;
        std::cerr << "Should NOT be isomorphic" << std::endl;
        assert(false);
      }
    }
  } 
}


int main(int, char**)
{
  testTree(tree1); 
  testTree(tree2); 
  testTree(treeWithDistances1);
  testDistances(true);
  testDistances(false);
  testMADRooting();
  testIsomorphism();
  std::cout << "Test PLLUnrootedTree ok!" << std::endl;
  return 0;
}







