#pragma once

#include <IO/FamiliesFileParser.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <trees/PLLUnrootedTree.hpp>

class PerCoreGeneTrees {
public:
  struct GeneTree {
    
    std::string name;
    GeneSpeciesMapping mapping;
    PLLUnrootedTree *geneTree;
    bool ownTree; // If true, I am responsible for destroying the tree
    ~GeneTree() {
      if (ownTree) {
        delete geneTree;
      }
    }
  };

  PerCoreGeneTrees(const Families &families);
  PerCoreGeneTrees(const GeneSpeciesMapping &mapping, PLLUnrootedTree &geneTree);
  
  PerCoreGeneTrees(const PerCoreGeneTrees &) = delete;
  PerCoreGeneTrees & operator = (const PerCoreGeneTrees &) = delete;
  PerCoreGeneTrees(PerCoreGeneTrees &&) = delete;
  PerCoreGeneTrees & operator = (PerCoreGeneTrees &&) = delete;
  
  
  std::vector<GeneTree> &getTrees() {return _geneTrees;}
  bool checkMappings(const std::string &speciesTreeFile);
  double getLoadBalanceRatio() const {return _loadBalanceRatio;}
private:
  std::vector<GeneTree> _geneTrees;
  std::vector<unsigned int> _treeSizes;
  double _loadBalanceRatio;
};
