#pragma once

#include <IO/FamiliesFileParser.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/LibpllEvaluation.hpp>


class PerCoreGeneTrees {
public:
  struct GeneTree {
    GeneSpeciesMapping mapping;
    pll_utree_t *tree;
  };

  PerCoreGeneTrees(const vector<FamiliesFileParser::FamilyInfo> &families);
  
private:
  vector<GeneTree> _geneTrees;

};
