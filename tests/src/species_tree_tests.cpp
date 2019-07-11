#include <trees/SpeciesTree.hpp>
#include <cassert>

void checkMove(SpeciesTree &speciesTree, unsigned int direction) 
{
  std::string initialStr = speciesTree.toString();
  auto initialTaxa = speciesTree.getTaxaNumber();
  
  SpeciesTreeOperator::changeRoot(speciesTree, direction);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  assert(initialStr == speciesTree.toString());
}


void testRootMoves() 
{
  std::string initialTreeStr = "((A,B),(C,D));";
  SpeciesTree speciesTree(initialTreeStr, false);
  for (unsigned int i = 0; i < 4; ++i) {
    checkMove(speciesTree, i);
  }
  std::cout << "Test species tree ok!" << std::endl;
}

void testBuildRandomTree()
{
  std::unordered_set<std::string> labels;
  for (unsigned int i = 0; i < 20; ++i) {
    labels.insert(std::string("S" + std::to_string(i)));
  }
  SpeciesTree speciesTree(labels);
  std::cout << speciesTree << std::endl;
}

int main(int argc, char** argv)
{
  testRootMoves();
  testBuildRandomTree();

  return 0;
}
