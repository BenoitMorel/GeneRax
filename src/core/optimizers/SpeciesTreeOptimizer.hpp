#pragma once

#include <trees/SpeciesTree.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <string>
#include <maths/Parameters.hpp>
#include <maths/AverageStream.hpp>
#include <util/enums.hpp>
#include <IO/Families.hpp>
#include <memory>
#include <maths/ModelParameters.hpp>
#include <trees/Clade.hpp>



struct EvaluatedMove {
  unsigned int prune;
  unsigned int regraft;
  double ll;
};

struct DistanceInfo {
  DistanceMatrix distanceMatrix;
  std::vector<std::string> speciesIdToSpeciesString;
  StringToUint speciesStringToSpeciesId;
};

struct MovesBlackList;

class SpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  SpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &initialFamilies, 
      const RecModelInfo &recModelInfo,
      const Parameters &startingRates,
      bool userDTLRates,
      const std::string &outputDir,
      bool constrainSearch);
  
  // forbid copy
  SpeciesTreeOptimizer(const SpeciesTreeOptimizer &) = delete;
  SpeciesTreeOptimizer & operator = (const SpeciesTreeOptimizer &) = delete;
  SpeciesTreeOptimizer(SpeciesTreeOptimizer &&) = delete;
  SpeciesTreeOptimizer & operator = (SpeciesTreeOptimizer &&) = delete;
  virtual ~SpeciesTreeOptimizer(); 
   
  enum OptimizationCriteria {
    ReconciliationLikelihood = 0,
    SupportedClades
  };

  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);

  void optimize(SpeciesSearchStrategy strategy,
      unsigned int sprRadius,
      OptimizationCriteria criteria = ReconciliationLikelihood);
  
  double optimizeDTLRates();
  
  double rootSearch(unsigned int maxDepth = 1000000);

  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
  const SpeciesTree &getSpeciesTree() const {return *_speciesTree;} 

  double getReconciliationLikelihood() const {return _bestRecLL;}

  double computeRecLikelihood();

private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  std::unique_ptr<PerCoreGeneTrees> _geneTrees;
  PerCoreEvaluations _evaluations; 
  std::vector<pll_unode_t*> _previousGeneRoots;
  Families _initialFamilies;
  std::string _outputDir;
  double _lastRecLL;
  double _bestRecLL;
  bool _firstOptimizeRatesCall;
  bool _userDTLRates;
  ModelParameters _modelRates;
  CladeSet _geneClades;
  bool _constrainSearch;
  unsigned int _okForClades;
  unsigned int _koForClades;
  AverageStream _averageGeneRootDiff;
  bool _hardToFindBetter;
  OptimizationCriteria _optimizationCriteria;
  DistanceInfo _distanceInfo;
private:
  void _computeAllGeneClades();
  unsigned int _unsupportedCladesNumber();
  void _computeDistanceInfo();
  ModelParameters computeOptimizedRates(); 
  void updateEvaluations();
  void rootSearchAux(SpeciesTree &speciesTree, 
      PerCoreGeneTrees &geneTrees, 
      RecModel model, 
      std::vector<unsigned int> &movesHistory, 
      std::vector<unsigned int> &bestMovesHistory, 
      double &bestLL, 
      unsigned int &visits,
      unsigned int maxDepth);
  bool testPruning(unsigned int prune,
    unsigned int regraft);
  void newBestTreeCallback();
  void beforeTestCallback();
  void rollbackCallback();
  std::string getSpeciesTreePath(const std::string &speciesId);
  void setGeneTreesFromFamilies(const Families &families);
  void reGenerateEvaluations();
  void optimizeGeneRoots();
  double transferSearch();
  double transferRound(MovesBlackList &blacklist, 
      bool &maxImprovementsReached);
  double sprSearch(unsigned int radius);
  double fastSPRRound(unsigned int radius);
  double veryLocalSearch(unsigned int spid);
  void setOptimizationCriteria(OptimizationCriteria criteria) {
    _optimizationCriteria = criteria;
  }
};
