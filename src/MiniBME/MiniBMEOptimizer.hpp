#pragma once

#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Families.hpp>
#include <memory>
#include <vector>
#include <NJ/MiniBME.hpp>
#include <NJ/MiniBMEPruned.hpp>



class USearchMiniBMEEvaluator {
public:
  USearchMiniBMEEvaluator(PLLUnrootedTree &speciesTree,
    const Families &families,
    bool missingData):
      _lastScore(-99999) {
      if (missingData) {
        _miniBME.reset(new MiniBMEPruned(speciesTree, families)); 
      } else {
        _miniBME.reset(new MiniBME(speciesTree, families)); 
      }
    }

  virtual ~USearchMiniBMEEvaluator() {}
  virtual double eval(PLLUnrootedTree &tree);
  bool computeAndApplyBestSPR(PLLUnrootedTree &tree,
    unsigned int maxRadiusWithoutImprovement);
private:
  std::unique_ptr<BMEEvaluator> _miniBME;
  double _lastScore;
};


class MiniBMEOptimizer {
public:
  MiniBMEOptimizer(const std::string speciesTreeFile, 
      const Families &families,
      bool missingData,
      const std::string &outputDir);

  void optimize();
private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  std::string _outputDir;
  bool _missingData;
  Families _families;
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



