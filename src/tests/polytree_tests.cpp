#include <iostream>

#include <trees/PLLUnrootedTree.hpp>
#include <trees/PolyTree.hpp>

void test_1()
{
  std::string poly1 = "(A_1:1.0, (A_2:1.0, B_2:1.0):0.000001):1.0";
  std::string poly2 = "((A_3:1.0, B_3:1.0):0.000001, (C_3:1.0,D_3:1.0):0.000001):1.0";
  std::string poly3 = "((A_5:1.0, B_5:1.0):0.000001, (C_5:1.0,D_5:1.0):0.000001):1.0";
  std::string geneTreeStr = "(" + poly1 + "," + poly2 + "," + poly3 + ");";
  PLLUnrootedTree geneTree(geneTreeStr, false); 
  PolyTree polyTree(geneTree);
  std::cout << "Gene tree no   polytomy: " << std::endl;
  std::cout << geneTree.getNewickString() << std::endl;
  std::cout << "Gene tree with polytomy: " << std::endl;
  std::cout << polyTree.getNewickString() << std::endl;
}

int main(int, char**)
{
  test_1();
  return 0;
}



