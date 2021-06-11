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
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <search/SpeciesRootSearch.hpp>
#include <trees/Clade.hpp>
#include <util/Constants.hpp>
#include <util/types.hpp>


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

struct SpeciesTreeSearchParams {
  SpeciesTreeSearchParams():
    sprRadius(DEFAULT_SPECIES_SPR_RADIUS),
    rootSmallRadius(DEFAULT_SPECIES_SMALL_ROOT_RADIUS),
    rootBigRadius(DEFAULT_SPECIES_BIG_ROOT_RADIUS)
  {}
  unsigned int sprRadius;
  unsigned int rootSmallRadius;
  unsigned int rootBigRadius;
};

struct MovesBlackList;


class SpeciesTreeLikelihoodEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  SpeciesTreeLikelihoodEvaluator(PerCoreEvaluations &evaluations,
      bool rootedGeneTrees):
    _evaluations(evaluations),
    _rootedGeneTrees(rootedGeneTrees)
  {}
  virtual ~SpeciesTreeLikelihoodEvaluator() {}
  virtual double computeLikelihood();
  virtual void forceGeneRootOptimization();
  virtual void pushRollback();
  virtual void popAndApplyRollback();
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL);
private:
  PerCoreEvaluations &_evaluations;
  std::stack<std::vector<pll_unode_t*> > _previousGeneRoots;
  bool _rootedGeneTrees;
};


class SpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  SpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &initialFamilies, 
      const RecModelInfo &recModelInfo,
      const Parameters &startingRates,
      bool userDTLRates,
      const std::string &outputDir,
      const SpeciesTreeSearchParams &searchParams);
  
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
      OptimizationCriteria criteria = ReconciliationLikelihood);
  
  double optimizeDTLRates(bool thorough = false);
  
  double rootSearch(unsigned int maxDepth,
      bool optimizeParams,
      bool outputConsel);

  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
  const SpeciesTree &getSpeciesTree() const {return *_speciesTree;} 

  double getReconciliationLikelihood() const {return _bestRecLL;}

  double computeRecLikelihood();

  void addPerFamilyLikelihoods(const std::string &newick,
      TreePerFamLLVec &treePerFamLLVec);

  void savePerFamilyLikelihoods(const TreePerFamLLVec &treePerFamLLVec,
      const std::string &treesOutput,
      const std::string &llOutput);

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
  SpeciesTreeSearchParams _searchParams;
  unsigned int _okForClades;
  unsigned int _koForClades;
  AverageStream _averageGeneRootDiff;
  bool _hardToFindBetter;
  OptimizationCriteria _optimizationCriteria;
private:
  void _computeAllGeneClades();
  unsigned int _unsupportedCladesNumber();
  ModelParameters computeOptimizedRates(bool thorough); 
  void updateEvaluations();
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
  std::vector<double> _getSupport();
  
};
