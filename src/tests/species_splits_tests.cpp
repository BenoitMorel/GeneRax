#include <iostream>

#include <ccp/SpeciesSplits.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLUnrootedTree.hpp>


int main()
{
  std::string speciesTreeStr = "((A,B),(C,D),(E,F));";
  PLLUnrootedTree speciesTree(speciesTreeStr, false);
  std::vector<std::string> geneTreeStrings;
  geneTreeStrings.push_back("((A_1,B_1),((C_1,D_1),(C_2,D_2)), (E_1, F_1));");
  geneTreeStrings.push_back("((A_1,B_1),C_1, (E_1, F_1));");
  geneTreeStrings.push_back("((A_1,C_1),B_1, F_1);");

  SpeciesSplits speciesSplits(speciesTree.getLabels());

  for (auto &geneTreeStr: geneTreeStrings) {
    PLLUnrootedTree geneTree(geneTreeStr, false);
    GeneSpeciesMapping mapping;
    mapping.fillFromGeneLabels(geneTree.getLabels());
    speciesSplits.addGeneTree(geneTree, mapping);
    //std::cout << "Total number of distinct splits: " << speciesSplits.distinctSplitsNumber() << std::endl;
    //std::cout << "Total number of splits: " << speciesSplits.nonDistinctSplitsNumber() << std::endl;
  }
  assert(speciesSplits.distinctSplitsNumber() == 9);
  assert(speciesSplits.nonDistinctSplitsNumber() == 12);
  std::cout << "Test SpeciesSplits ok!" << std::endl;
  return 0;
}
