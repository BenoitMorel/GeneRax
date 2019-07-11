#include <trees/SpeciesTree.hpp>
#include <cassert>

void checkRootMove(SpeciesTree &speciesTree, unsigned int direction) 
{
  std::string initialStr = speciesTree.toString();
  auto initialTaxa = speciesTree.getTaxaNumber();
  
  SpeciesTreeOperator::changeRoot(speciesTree, direction);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  assert(initialStr == speciesTree.toString());
}

void checkSPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft) 
{
  std::string initialStr = speciesTree.toString();
  auto initialTaxa = speciesTree.getTaxaNumber();
  unsigned int rollback = SpeciesTreeOperator::applySPRMove(speciesTree, prune, regraft);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  SpeciesTreeOperator::reverseSPRMove(speciesTree, prune, rollback);
  assert(initialTaxa == speciesTree.getTaxaNumber());
  assert(initialStr == speciesTree.toString());
}

void testRootMoves() 
{
  std::string initialTreeStr = "((A,B),(C,D));";
  SpeciesTree speciesTree(initialTreeStr, false);
  for (unsigned int i = 0; i < 4; ++i) {
    checkRootMove(speciesTree, i);
  }
  std::cout << "Test species tree ok!" << std::endl;
}

void testSPRMoves()
{
  std::string initialTreeStr = "((A, (B, C)),((D, E), (F, G)));";
  SpeciesTree speciesTree(initialTreeStr, false);
  unsigned int radius = 10;
  for (unsigned int i = 0; i < speciesTree.getMaxNodeIndex(); ++i) {
    std::vector<unsigned int> regrafts;
    auto pruneNode = speciesTree.getNode(i);
    if (pruneNode == speciesTree.getRoot()) {
      continue;
    }
    SpeciesTreeOperator::getPossibleRegrafts(speciesTree, pruneNode->node_index, radius, regrafts);
    for (auto regraft: regrafts) {
      checkSPRMove(speciesTree, pruneNode->node_index, regraft); 
    }
  }
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
  testSPRMoves();
  return 0;
}

