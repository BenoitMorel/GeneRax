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
#include <search/SpeciesSearchCommon.hpp>

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
  SpeciesTreeLikelihoodEvaluator() {}
  void init(PerCoreEvaluations &evaluations,
      PerCoreGeneTrees &geneTrees,
      ModelParameters &modelRates,
      bool rootedGeneTrees,
      bool pruneSpeciesTree,
      bool userDTLRates) 
  {
    _evaluations = &evaluations;
    _geneTrees = &geneTrees;
    _modelRates = &modelRates;
    _rootedGeneTrees = rootedGeneTrees;
    _pruneSpeciesTree = pruneSpeciesTree;
    _userDTLRates = userDTLRates;
  }
  virtual ~SpeciesTreeLikelihoodEvaluator() {}
  virtual double computeLikelihood();
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const;
  virtual bool isDated() const {return _modelRates->info.isDated();}
  virtual double optimizeModelRates(bool thorough = false);
  virtual void pushRollback();
  virtual void popAndApplyRollback();
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL);
  virtual void getTransferInformation(SpeciesTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents);
  virtual bool pruneSpeciesTree() const {return _pruneSpeciesTree;}
private:
  PerCoreGeneTrees *_geneTrees;
  PerCoreEvaluations *_evaluations;
  ModelParameters *_modelRates;
  std::stack<std::vector<corax_unode_t*> > _previousGeneRoots;
  bool _rootedGeneTrees;
  bool _pruneSpeciesTree;
  bool _userDTLRates;
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
   
  virtual void onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);

  void optimize(SpeciesSearchStrategy strategy);
  
  double rootSearch(unsigned int maxDepth,
      bool outputConsel);

  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
  const SpeciesTree &getSpeciesTree() const {return *_speciesTree;} 

  double getReconciliationLikelihood() const {return _searchState.bestLL;}

  double computeRecLikelihood();

  void savePerFamilyLikelihoods(const TreePerFamLLVec &treePerFamLLVec,
      const std::string &treesOutput,
      const std::string &llOutput);

private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  std::unique_ptr<PerCoreGeneTrees> _geneTrees;
  PerCoreEvaluations _evaluations; 
  SpeciesTreeLikelihoodEvaluator _evaluator;
  std::vector<corax_unode_t*> _previousGeneRoots;
  Families _initialFamilies;
  std::string _outputDir;
  bool _firstOptimizeRatesCall;
  bool _userDTLRates;
  ModelParameters _modelRates;
  CladeSet _geneClades;
  SpeciesTreeSearchParams _searchParams;
  unsigned int _okForClades;
  unsigned int _koForClades;
  SpeciesSearchState _searchState;
private:
  void _computeAllGeneClades();
  unsigned int _unsupportedCladesNumber();
  void updateEvaluations();
  std::string getSpeciesTreePath(const std::string &speciesId);
  void setGeneTreesFromFamilies(const Families &families);
  void reGenerateEvaluations();
  double transferSearch();
  double sprSearch(unsigned int radius);
  std::vector<double> _getSupport();
  
};
