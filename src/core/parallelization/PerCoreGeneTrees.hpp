#pragma once

#include <IO/FamiliesFileParser.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <trees/PLLUnrootedTree.hpp>

/**
 * Holds the gene trees (with there mappings to the species tree)
 * allocated to the current parallel core.
 *
 * The sets of PerCoreGeneTrees over all parallel cores are a 
 * partition of all gene trees.
 */
class PerCoreGeneTrees {
public:
  struct GeneTree {
    
    std::string name;
    unsigned int familyIndex;
    GeneSpeciesMapping mapping;
    PLLUnrootedTree *geneTree;
    bool ownTree; // If true, I am responsible for destroying the tree
    ~GeneTree() {
      if (ownTree) {
        delete geneTree;
      }
    }
  };

  /**
   *  Parse and allocate to the current core the gene trees 
   *  from the description of the gene families.
   *  @param families families description
   */
  PerCoreGeneTrees(const Families &families, 
      bool acceptMultipleTrees = false);
  /**
   * Create an instance with a unique gene tree, without 
   * accouting for parallelization.
   *
   * This is useful for applying some methods, whose interface
   * require a PerCoreGeneTrees object, to a unique gene tree
   * without parallelization.
   *  @param mapping The gene to species mapping
   *  @param geneTree The gene tree
   */
  PerCoreGeneTrees(const GeneSpeciesMapping &mapping, PLLUnrootedTree &geneTree);
  
  PerCoreGeneTrees(const PerCoreGeneTrees &) = delete;
  PerCoreGeneTrees & operator = (const PerCoreGeneTrees &) = delete;
  PerCoreGeneTrees(PerCoreGeneTrees &&) = delete;
  PerCoreGeneTrees & operator = (PerCoreGeneTrees &&) = delete;
  
  /**
   *  @return Trees allocated to the current core
   */
  std::vector<GeneTree> &getTrees() {return _geneTrees;}
  const std::vector<GeneTree> &getTrees() const {return _geneTrees;}

  /**
   *  @param speciesTreeFile path to the species tree file
   *  @return true if the mappings are valid.
   */
  bool checkMappings(const std::string &speciesTreeFile);

  static void getPerCoreFamilies(const Families &allFamilies,
      Families &perCoreFamilies);

private:
  std::vector<GeneTree> _geneTrees;
  std::vector<unsigned int> _treeSizes;

};
