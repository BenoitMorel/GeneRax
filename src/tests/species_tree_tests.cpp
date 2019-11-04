#include <trees/SpeciesTree.hpp>
#include <cassert>

void checkRootMove(SpeciesTree &speciesTree, unsigned int direction) 
{
  std::string initialStr = speciesTree.toString();
  auto initialTaxa = speciesTree.getTree().getLeavesNumber();
  
  SpeciesTreeOperator::changeRoot(speciesTree, direction);
  assert(initialTaxa == speciesTree.getTree().getLeavesNumber());
  SpeciesTreeOperator::revertChangeRoot(speciesTree, direction);
  assert(initialTaxa == speciesTree.getTree().getLeavesNumber());
  assert(initialStr == speciesTree.toString());
}

void checkSPRMove(SpeciesTree &speciesTree, unsigned int prune, unsigned int regraft) 
{
  std::string initialStr = speciesTree.toString();
  auto initialTaxa = speciesTree.getTree().getLeavesNumber();
  unsigned int rollback = SpeciesTreeOperator::applySPRMove(speciesTree, prune, regraft);
  assert(initialTaxa == speciesTree.getTree().getLeavesNumber());
  SpeciesTreeOperator::reverseSPRMove(speciesTree, prune, rollback);
  assert(initialTaxa == speciesTree.getTree().getLeavesNumber());
  assert(initialStr == speciesTree.toString());
}

void testRootMoves() 
{
  std::string initialTreeStr = "((A,B),(C,D));";
  SpeciesTree speciesTree(initialTreeStr, false);
  for (unsigned int i = 0; i < 4; ++i) {
    checkRootMove(speciesTree, i);
  }
}

void testSPRMoves()
{
  std::string initialTreeStr = "((A, (B, C)),((D, E), (F, G)));";
  SpeciesTree speciesTree(initialTreeStr, false);
  unsigned int radius = 10;
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(speciesTree, prunes);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      checkSPRMove(speciesTree, prune, regraft); 
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
}

int main(int, char**)
{
  testRootMoves();
  testBuildRandomTree();
  testSPRMoves();
  std::cout << "Test species tree ok!" << std::endl;
  return 0;
}

