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
      const RecModelInfo &recModelInfo,
      const Parameters &startingRates,
      bool userDTLRates,
      const std::string &outputDir,
      const std::string &execPath,
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
  
  double rootExhaustiveSearch();
  double transferSearch();
  double fastTransfersRound(MovesBlackList &blacklist);
  double reconciliationSearch();
  double reconciliationRound();
  double sprSearch(unsigned int radius);
  double fastSPRRound(unsigned int radius);
  double optimizeDTLRates();
 

  double computeLikelihood();
  // return the path to the saved species tree
  std::string saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
  const SpeciesTree &getSpeciesTree() const {return *_speciesTree;} 

  double getReconciliationLikelihood() const {return _bestRecLL;}

  double computeRecLikelihood();


  void likelihoodsSnapshot();


private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  std::unique_ptr<PerCoreGeneTrees> _geneTrees;
  PerCoreEvaluations _evaluations; 
  std::vector<pll_unode_t*> _previousGeneRoots;
  Families _initialFamilies;
  std::string _outputDir;
  std::string _execPath;
  double _lastRecLL;
  double _bestRecLL;
  SpeciesSearchStats _stats;
  bool _firstOptimizeRatesCall;
  bool _userDTLRates;
  ModelParameters _modelRates;
  CladeSet _geneClades;
  bool _constrainSearch;
  unsigned int _okForClades;
  unsigned int _koForClades;
  AverageStream _averageGeneRootDiff;
private:
  void _computeAllGeneClades();
  unsigned int _unsupportedCladesNumber();
  ModelParameters computeOptimizedRates(); 
  void updateEvaluations();
  void rootExhaustiveSearchAux(SpeciesTree &speciesTree, 
      PerCoreGeneTrees &geneTrees, 
      RecModel model, 
      std::vector<unsigned int> &movesHistory, 
      std::vector<unsigned int> &bestMovesHistory, 
      double &bestLL, 
      unsigned int &visits);
  bool testPruning(unsigned int prune,
    unsigned int regraft);
  void newBestTreeCallback();
  void beforeTestCallback();
  void rollbackCallback();
  std::string getSpeciesTreePath(const std::string &speciesId);
  void setGeneTreesFromFamilies(const Families &families);
  void reGenerateEvaluations();
  void optimizeGeneRoots();
};
