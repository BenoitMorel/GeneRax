#pragma once

#include <trees/SpeciesTree.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <string>
#include <maths/Parameters.hpp>
#include <util/enums.hpp>
#include <families/Families.hpp>
#include <memory>
  
struct EvaluatedMove {
  unsigned int prune;
  unsigned int regraft;
  double ll;
};

class SpeciesTreeOptimizer: public SpeciesTree::Listener {
public:
  SpeciesTreeOptimizer(const std::string speciesTreeFile, 
      const Families &initialFamilies, 
      RecModel model,
      double supportThreshold,
      const std::string &outputDir,
      const std::string &execPath);
  
  // forbid copy
  SpeciesTreeOptimizer(const SpeciesTreeOptimizer &) = delete;
  SpeciesTreeOptimizer & operator = (const SpeciesTreeOptimizer &) = delete;
  SpeciesTreeOptimizer(SpeciesTreeOptimizer &&) = delete;
  SpeciesTreeOptimizer & operator = (SpeciesTreeOptimizer &&) = delete;
  virtual ~SpeciesTreeOptimizer(); 
    
  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);

  void rootExhaustiveSearch(bool doOptimizeGeneTrees);
  double fastSPRRound(unsigned int radius);
  double slowSPRRound(unsigned int radius, double bestLL);
  double sprSearch(unsigned int radius, bool doOptimizeGeneTrees);
  void optimizeDTLRates();
 
  Parameters computeOptimizedRates() const; 

  double optimizeGeneTrees(unsigned int radius);
  void revertGeneTreeOptimization();

  double computeLikelihood(unsigned int geneSPRRadius);
  void saveCurrentSpeciesTreeId(std::string str = "inferred_species_tree.newick", bool masterRankOnly = true);
  void saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly = true);
  const SpeciesTree &getSpeciesTree() const {return *_speciesTree;} 

  double getReconciliationLikelihood() const {return _bestRecLL;}
  double getLibpllLikeliohood() const {return _bestLibpllLL;}
  double getJointLikelihood() const {return getReconciliationLikelihood() + getLibpllLikeliohood();}

  double computeReconciliationLikelihood();
private:
  std::unique_ptr<SpeciesTree> _speciesTree;
  std::unique_ptr<PerCoreGeneTrees> _geneTrees;
  std::vector<std::shared_ptr<ReconciliationEvaluation> > _evaluations; // one per geneTree
  Families _initialFamilies;
  Families _currentFamilies;
  RecModel _model;
  std::string _outputDir;
  std::string _execPath;
  unsigned int _geneTreeIteration;
  Parameters _hack;
  double _supportThreshold;
  double _lastRecLL;
  double _lastLibpllLL;
  double _bestRecLL;
  double _bestLibpllLL;
private:
  void updateEvaluations();
  void rootExhaustiveSearchAux(SpeciesTree &speciesTree, 
      PerCoreGeneTrees &geneTrees, 
      RecModel model, 
      bool doOptimizeGeneTrees, 
      std::vector<unsigned int> &movesHistory, 
      std::vector<unsigned int> &bestMovesHistory, 
      double &bestLL, 
      unsigned int &visits);
  void newBestTreeCallback();
  std::string getSpeciesTreePath(const std::string &speciesId);
  std::vector<EvaluatedMove> getSortedCandidateMoves(unsigned int speciesRadius);
  void setGeneTreesFromFamilies(const Families &families);
  void reGenerateEvaluations();
};
