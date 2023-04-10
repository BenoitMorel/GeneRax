#include <iostream>
#include <trees/PLLRootedTree.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <trees/PolytomySolver.hpp>

void test_simple_interface_1()
{
  PLLRootedTree speciesTree("((A,B)E,(C,D)F)G;", false);
  std::map<std::string, unsigned int> speciesLabelsToSolve;
  speciesLabelsToSolve.insert({std::string("A"), 4});
  speciesLabelsToSolve.insert({std::string("B"), 2});
  speciesLabelsToSolve.insert({std::string("C"), 1});
  speciesLabelsToSolve.insert({std::string("E"), 1});
  PolytomySolver::solveSimpleInterface(speciesTree,
      speciesLabelsToSolve);
}

void test_simple_interface_2()
{
  PLLRootedTree speciesTree("((A,B)E,(C,D)F)G;", false);
  std::map<std::string, unsigned int> speciesLabelsToSolve;
  speciesLabelsToSolve.insert({std::string("A"), 2});
  speciesLabelsToSolve.insert({std::string("B"), 2});
  speciesLabelsToSolve.insert({std::string("C"), 2});
  speciesLabelsToSolve.insert({std::string("D"), 2});
  PolytomySolver::solveSimpleInterface(speciesTree,
      speciesLabelsToSolve);
}

void test_simple_interface_3()
{
  PLLRootedTree speciesTree("((A,B)E,(C,D)F)G;", false);
  std::map<std::string, unsigned int> speciesLabelsToSolve;
  speciesLabelsToSolve.insert({std::string("E"), 1});
  speciesLabelsToSolve.insert({std::string("F"), 1});
  speciesLabelsToSolve.insert({std::string("G"), 1});
  PolytomySolver::solveSimpleInterface(speciesTree,
      speciesLabelsToSolve);
}

void test_simple_interface_4()
{
  PLLRootedTree speciesTree("((A,B)E,(C,D)F)G;", false);
  std::map<std::string, unsigned int> speciesLabelsToSolve;
  speciesLabelsToSolve.insert({std::string("A"), 2});
  speciesLabelsToSolve.insert({std::string("B"), 2});
  speciesLabelsToSolve.insert({std::string("C"), 1});
  speciesLabelsToSolve.insert({std::string("D"), 1});
  speciesLabelsToSolve.insert({std::string("F"), 1});
  PolytomySolver::solveSimpleInterface(speciesTree,
      speciesLabelsToSolve);
}

void test_1()
{
  /*
  PLLRootedTree speciesTree("((A,B),(C,D));", false);
  std::string poly1 = "(A_1:1.0, (A_2:1.0, B_2:1.0):0.000001):1.0";
  std::string poly2 = "((A_3:1.0, B_3:1.0):0.000001, (C_3:1.0,D_3:1.0):0.000001):1.0";
  std::string poly3 = "((A_5:1.0, B_5:1.0):0.000001, (C_5:1.0,D_5:1.0):0.000001):1.0";
  std::string geneTreeStr = "(" + poly1 + "," + poly2 + "," + poly3 + ");";
  PLLUnrootedTree geneTree(geneTreeStr, false);
  
  auto speciesLabelToId = speciesTree.getLabelToIntMap();
  std::vector<unsigned int> geneToSpecies(
      geneTree.getDirectedNodeNumber(), 0);
  for (auto leaf: geneTree.getLeaves()) {
    std::string geneString = leaf->label;
    std::string speciesString;
    speciesString += geneString[0];
    geneToSpecies[leaf->node_index] = speciesLabelToId[speciesString];
  }
  PolytomySolver::solve(speciesTree, geneTree, geneToSpecies);
*/
}

int main(int, char**)
{
  test_simple_interface_1();
  test_simple_interface_2();
  test_simple_interface_3();
  test_simple_interface_4();
  return 0;
}


