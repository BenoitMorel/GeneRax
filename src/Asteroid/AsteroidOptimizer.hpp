#pragma once

#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include <IO/Families.hpp>
#include <memory>
#include <vector>
#include <DistanceMethods/Asteroid.hpp>
#include <DistanceMethods/Astrid.hpp>



class USearchAsteroidEvaluator {
public:
  USearchAsteroidEvaluator(PLLUnrootedTree &speciesTree,
    const Families &families,
    double minbl,
    bool missingData):
      _lastScore(-99999) {
      if (missingData) {
        _miniBME.reset(new Asteroid(speciesTree, families, minbl)); 
      } else {
        _miniBME.reset(new Astrid(speciesTree, families, minbl)); 
      }
    }

  virtual ~USearchAsteroidEvaluator() {}
  virtual double eval(PLLUnrootedTree &tree);
  bool computeAndApplyBestSPR(PLLUnrootedTree &tree,
    unsigned int maxRadiusWithoutImprovement);
private:
  std::unique_ptr<BMEEvaluator> _miniBME;
  double _lastScore;
};


class AsteroidOptimizer {
public:
  AsteroidOptimizer(const std::string speciesTreeFile, 
      const Families &families,
      double minbl,
      bool missingData,
      const std::string &outputDir);

  void optimize();
private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  std::string _outputDir;
  double _minbl;
  bool _missingData;
  Families _families;
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
};



