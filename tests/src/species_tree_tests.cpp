#include <trees/SpeciesTree.hpp>
#include <cassert>

void checkMove(SpeciesTree &speciesTree, bool left1, bool left2) 
{
  std::string initialStr = speciesTree.toString();
  auto initialTaxa = speciesTree.getTaxaNumber();
  
  SpeciesTreeOperator::changeRoot(speciesTree, left1, left2);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  SpeciesTreeOperator::changeRoot(speciesTree, !left1, !left2);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  assert(initialStr == speciesTree.toString());
}

int main(int argc, char** argv)
{
  std::string initialTreeStr = "((A,B),(C,D));";
  SpeciesTree speciesTree(initialTreeStr, false);

  checkMove(speciesTree, true, true);
  checkMove(speciesTree, true, false);
  checkMove(speciesTree, false, true);
  checkMove(speciesTree, false, false);
  std::cout << "Test species tree ok!" << std::endl;
  return 0;
}

