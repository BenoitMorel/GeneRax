#pragma once

#include <trees/SpeciesTree.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <string>
#include <maths/Parameters.hpp>
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

struct MovesBlackList;

struct SpeciesSearchStats {
  unsigned int exactLikelihoodCalls;
  unsigned int approxLikelihoodCalls;
  unsigned int testedTrees;
  unsigned int testedTransfers;
  unsigned int acceptedTrees;
  unsigned int acceptedTransfers;
  SpeciesSearchStats() { reset(); }

  friend std::ostream& operator<<(std::ostream& os , const SpeciesSearchStats &) {
    /*
    os << "Tested trees: " << stats.testedTrees << std::endl;
    os << "Accepted trees: " << stats.acceptedTrees << std::endl;
    os << "Tested transfers: " << stats.testedTransfers << std::endl;
    os << "Accepted transfers trees: " << stats.acceptedTransfers << std::endl;
    */
    //os << "Approx likelihood calls: " << stats.approxLikelihoodCalls << std::endl;
    //os << "Exact likelihood calls: " << stats.exactLikelihoodCalls << std::endl;
    return os;
  }
  void reset() {
    exactLikelihoodCalls = 0;
    approxLikelihoodCalls = 0;
    testedTrees = 0;
    testedTransfers = 0;
    acceptedTrees = 0;
    acceptedTransfers = 0;
  }
};

class SpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  SpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &initialFamilies, 
      RecModel model,
      const Parameters &startingRates,
      bool perFamilyRates,
      bool userDTLRates,
      double minGeneBranchLength,
      bool pruneSpeciesTree,
      double supportThreshold,
      const std::string &outputDir,
      const std::string &execPath,
      bool fractionMissing,
      bool constrainSearch);
  
  // forbid copy
  SpeciesTreeOptimizer(const SpeciesTreeOptimizer &) = delete;
  SpeciesTreeOptimizer & operator = (const SpeciesTreeOptimizer &) = delete;
  SpeciesTreeOptimizer(SpeciesTreeOptimizer &&) = delete;
  SpeciesTreeOptimizer & operator = (SpeciesTreeOptimizer &&) = delete;
  virtual ~SpeciesTreeOptimizer(); 
   


  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);

  void optimize(SpeciesSearchStrategy strategy,
      unsigned int sprRadius);
  
  double rootExhaustiveSearch(bool doOptimizeGeneTrees);
  double fastTransfersRound(MovesBlackList &blacklist);
  double fastSPRRound(unsigned int radius);
  double slowSPRRound(unsigned int radius, double bestLL);
  double sprSearch(unsigned int radius, bool doOptimizeGeneTrees = false);
  double transferSearch();
  double optimizeDTLRates();
 

  double optimizeGeneTrees(unsigned int radius);
  void revertGeneTreeOptimization();

  double computeLikelihood(unsigned int geneSPRRadius);
  // return the path to the saved species tree
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
  const SpeciesTree &getSpeciesTree() const {return *_speciesTree;} 

  double getReconciliationLikelihood() const {return _bestRecLL;}
  double getLibpllLikeliohood() const {return _bestLibpllLL;}
  double getJointLikelihood() const {return getReconciliationLikelihood() + getLibpllLikeliohood();}

  double computeRecLikelihood();
  double computeApproxRecLikelihood();


  void likelihoodsSnapshot();


private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  std::unique_ptr<PerCoreGeneTrees> _geneTrees;
  PerCoreEvaluations _evaluations; 
  Families _initialFamilies;
  Families _currentFamilies;
  RecModel _recModel;
  std::string _outputDir;
  std::string _execPath;
  unsigned int _geneTreeIteration;
  double _supportThreshold;
  double _lastRecLL;
  double _lastLibpllLL;
  double _bestRecLL;
  double _bestLibpllLL;
  SpeciesSearchStats _stats;
  bool _firstOptimizeRatesCall;
  bool _userDTLRates;
  double _minGeneBranchLength;
  bool _pruneSpeciesTree;
  ModelParameters _modelRates;
  std::string _fractionMissingFile;
  CladeSet _geneClades;
  bool _constrainSearch;
  unsigned int _okForClades;
  unsigned int _koForClades;
private:
  void _computeAllGeneClades();
  unsigned int _unsupportedCladesNumber();
  ModelParameters computeOptimizedRates(); 
  void updateEvaluations();
  void rootExhaustiveSearchAux(SpeciesTree &speciesTree, 
      PerCoreGeneTrees &geneTrees, 
      RecModel model, 
      bool doOptimizeGeneTrees, 
      std::vector<unsigned int> &movesHistory, 
      std::vector<unsigned int> &bestMovesHistory, 
      double &bestLL, 
      unsigned int &visits);
  bool testPruning(unsigned int prune,
    unsigned int regraft,
    double refApproxLL, 
    unsigned int hash1);
  void newBestTreeCallback();
  std::string getSpeciesTreePath(const std::string &speciesId);
  std::vector<EvaluatedMove> getSortedCandidateMoves(unsigned int speciesRadius);
  void setGeneTreesFromFamilies(const Families &families);
  void reGenerateEvaluations();
};
