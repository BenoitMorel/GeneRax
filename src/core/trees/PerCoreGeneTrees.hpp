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
    bool ownTree; // If true, I am responsible for destroying the tree
    ~GeneTree() {
      if (ownTree) {
        pll_utree_destroy(tree, 0);}
    }
  };

  PerCoreGeneTrees(const Families &families);
  PerCoreGeneTrees(const GeneSpeciesMapping &mapping, pll_utree_t *tree);
  std::vector<GeneTree> &getTrees() {return _geneTrees;}
  bool checkMappings(const std::string &speciesTreeFile);
private:
  std::vector<GeneTree> _geneTrees;
  std::vector<unsigned int> _treeSizes;

};
