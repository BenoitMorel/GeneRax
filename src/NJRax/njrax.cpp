#include <string>
#include <iostream>

#include <DistanceMethods/MiniNJ.hpp>
#include <DistanceMethods/Cherry.hpp>
#include <DistanceMethods/CherryPro.hpp>
#include <trees/PLLRootedTree.hpp>
#include <IO/Families.hpp>
#include <util/enums.hpp>
#include <routines/SlavesMain.hpp>
#include <routines/Routines.hpp>

/**
 *  Hack for fix a link error
 */
int unused(int argc, char** argv)
{
  return static_scheduled_main(argc, argv, 0);
}

int main(int argc, char** argv)
{

  if (argc != 5) {
    std::cerr << "Error: syntax is ./njrax algorithm input_gene_trees mapping_file output_species_tree " << std::endl;
    std::cerr << "Possible algorithms are [MiniNJ, Cherry, NJst]" << std::endl;
    std::cerr << "The mapping file should have one of the two following syntaxes:" << std::endl;
   std::cerr << std::endl;
    std::cerr << "species1:gene1;gene2;gene3" << std::endl;
    std::cerr << "species2:gene4;gene6" << std::endl;
    std::cerr << std::endl;
    std::cerr << "or" << std::endl;
    std::cerr << std::endl;
    std::cerr << "gene1 species1" << std::endl;
    std::cerr << "gene2 species1" << std::endl;
    std::cerr << "gene3 species2" << std::endl;
    exit(1);
  }
  unsigned int a = 1;
  std::string algoStr(argv[a++]);
  std::string inputGeneTrees(argv[a++]);
  std::string inputMappingFile(argv[a++]);
  std::string outputSpeciesTree(argv[a++]);

  // let's hack to use GeneRax core interface
  Families families;
  FamilyInfo info;
  info.startingGeneTree = inputGeneTrees;
  info.mappingFile = inputMappingFile;
  families.push_back(info);
  SpeciesTreeAlgorithm speciesTreeAlgorithm = Enums::strToSpeciesTree(algoStr);
  std::string fakeOutputDir = "njrax_output";
  auto speciesTree(Routines::computeInitialSpeciesTree(families, 
        fakeOutputDir,
        speciesTreeAlgorithm));
  speciesTree->save(outputSpeciesTree);
  return 0;
}

