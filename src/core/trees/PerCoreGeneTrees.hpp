#pragma once

#include <IO/FamiliesFileParser.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/LibpllEvaluation.hpp>


class PerCoreGeneTrees {
public:
  struct GeneTree {
    std::string name;
    GeneSpeciesMapping mapping;
    std::string treeFile;
    pll_utree_t *tree;
    ~GeneTree() {pll_utree_destroy(tree, 0);}
  };

  PerCoreGeneTrees(const std::vector<FamiliesFileParser::FamilyInfo> &families);
  std::vector<GeneTree> &getTrees() {return _geneTrees;}
  bool checkMappings(const std::string &speciesTreeFile);
private:
  std::vector<GeneTree> _geneTrees;
  std::vector<unsigned int> _treeSizes;

};
