#include <iostream>

#include <ccp/SpeciesSplits.hpp>
#include <ccp/UnrootedSpeciesSplitScore.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include <trees/PLLUnrootedTree.hpp>


int main()
{
  std::string speciesTreeStr = "(((A,B),(C,D)),E);";
  PLLRootedTree speciesTree(speciesTreeStr, false);
  std::vector<std::string> geneTreeStrings;
  geneTreeStrings.push_back("((A_1,B_1),((C_1,D_1),(C_2,D_2)), E_1);");
  geneTreeStrings.push_back("((A_1,B_1),C_1, E_1);");
  geneTreeStrings.push_back("((A_1,C_1),B_1, E_1);");
  SpeciesSplits speciesSplits(speciesTree.getLabels(true), false);
  for (auto &geneTreeStr: geneTreeStrings) {
    PLLUnrootedTree geneTree(geneTreeStr, false);
    GeneSpeciesMapping mapping;
    mapping.fillFromGeneLabels(geneTree.getLabels());
    speciesSplits.addGeneTree(geneTree, mapping);
    //std::cout << "Total number of distinct splits: " << speciesSplits.distinctSplitsNumber() << std::endl;
    //std::cout << "Total number of splits: " << speciesSplits.nonDistinctSplitsNumber() << std::endl;
  }
 
  UnrootedSpeciesSplitScore score(speciesTree, speciesSplits); 
  std::cout << "Score: " << score.getScore() << std::endl;
  std::cout << "Test SplitScoreTest ok!" << std::endl;
  return 0;
}

