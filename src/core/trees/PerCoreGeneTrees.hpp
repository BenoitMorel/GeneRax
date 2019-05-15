#pragma once

#include <IO/FamiliesFileParser.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/LibpllEvaluation.hpp>


class PerCoreGeneTrees {
public:
  struct GeneTree {
    std::string name;
    GeneSpeciesMapping mapping;
    pll_utree_t *tree;
    ~GeneTree() {pll_utree_destroy(tree, 0);}
  };

  PerCoreGeneTrees(const std::vector<FamiliesFileParser::FamilyInfo> &families);
  const std::vector<GeneTree> &getTrees() {return _geneTrees;}
private:
  std::vector<GeneTree> _geneTrees;
  std::vector<unsigned int> _treeSizes;

};
