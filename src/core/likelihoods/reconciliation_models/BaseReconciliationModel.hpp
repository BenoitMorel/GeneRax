#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <util/Scenario.hpp>
#include <IO/Logger.hpp>
#include <util/enums.hpp>
#include <util/RecModelInfo.hpp>
#include <cmath>
#include <unordered_set>
#include <maths/ScaledValue.hpp>
#include <trees/PLLRootedTree.hpp>
#include <maths/Random.hpp>

#define IS_PROBA(x) (x == x) //(REAL(0.0) <= REAL(x) && REAL(x)  <= REAL(1.0)))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));

double log(ScaledValue v);
typedef std::vector< std::vector <double> > RatesVector;

inline double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}


/**
 *  Interface and common implementations for 
 *  all the reconciliation likelihood computation
 *  classes
 */
class ReconciliationModelInterface {
public:
  virtual ~ReconciliationModelInterface() {}
  
  /*
   * Set the per-species lineage rates
   */
  virtual void setRates(const RatesVector &rates) = 0;

  virtual void setAlpha(double) {}

  /**
   * (incrementally) compute and return the likelihood of the gene tree 
   */
  virtual double computeLogLikelihood() = 0;


  /**
   *  Sets the mode used to compute the incremental likelihood function
   */
  virtual void setPartialLikelihoodMode(PartialLikelihoodMode mode) = 0;
  
  /**
   * CLV invalidation for partial likelihood computation
   */
  virtual void invalidateAllSpeciesCLVs() = 0;
  /**
   *  Fill scenario with the maximum likelihood set of 
   *  events that would lead to the  current tree
   **/
  virtual bool inferMLScenario(Scenario &scenario, bool stochastic = false) = 0;

  /**
   *  Should be called after each change in the species tree topology
   */
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate) = 0;

  /**
   *  Set the path to the file containing information about
   *  fraction of missing genes
   */
  virtual void setFractionMissingGenes(
      const std::string &fractionMissingFile) = 0;
};

class BaseReconciliationModel: public ReconciliationModelInterface {
public:
  BaseReconciliationModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &recModelInfo); 

  virtual ~BaseReconciliationModel() {}
  
  // overload from parent
  virtual void invalidateAllSpeciesCLVs() {_allSpeciesNodesInvalid = true;}
  // overload from parent
  virtual void setPartialLikelihoodMode(PartialLikelihoodMode mode) {
    _likelihoodMode = mode;
  }
  // overload from parent
  virtual void setFractionMissingGenes(
      const std::string &fractionMissingFile);
  // overload from parent
  virtual void onSpeciesTreeChange(const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);

  corax_rnode_t *getSpeciesLeft(corax_rnode_t *node) {
    return _speciesLeft[node->node_index];
  }

  corax_rnode_t *getSpeciesRight(corax_rnode_t *node) {
    return _speciesRight[node->node_index];
  }

  corax_rnode_t *getSpeciesParent(corax_rnode_t *node) {
    return _speciesParent[node->node_index];
  }

  corax_rnode_t *getPrunedRoot() {return _prunedRoot;}
protected:

  /*
   * Init all structures that describe the species tree, in particular
   * the structure representing the pruned species tree in pruned mode
   * This function should be called once at the start
   */
  virtual void initSpeciesTree();
  
  /*
   *  Recompute probability values depending on the species tree only
   *  (not on the gene tree), such as the exctinction probability or
   *  the per-species event probabilities
   *  This function should be typically called after changing the DTL 
   *  rates or after updating the species tree.
   */
  virtual void recomputeSpeciesProbabilities() = 0;

  /**
   *  Callback that is always called at the start of computeLogLikelihood
   */
  virtual void beforeComputeLogLikelihood(); 
  
  /*
   * Fill nodes with the nodes of the species tree that belong to 
   * the "pruned species tree". nodesToAdd represents the species leaves
   * are covered by this gene family. Nodes are filled in post order fashion
   */
  bool fillPrunedNodesPostOrder(corax_rnode_t *node, 
    std::vector<corax_rnode_t *> &nodes, 
    std::unordered_set<corax_rnode_t *> *nodesToAdd = nullptr);  


  PLLRootedTree &getSpeciesTree() {return _speciesTree;}

  unsigned int getPrunedSpeciesNodeNumber() const {return _prunedSpeciesNodes.size();}
  unsigned int getAllSpeciesNodeNumber() const {return _allSpeciesNodes.size();}

  std::vector <corax_rnode_t *> &getAllSpeciesNodes() {return _allSpeciesNodes;}
  std::vector <corax_rnode_t *> &getPrunedSpeciesNodes() {return _prunedSpeciesNodes;}

  virtual size_t getSpeciesTreeHash() const;

protected:
  // description of the model
  RecModelInfo _info;
  // reference to the species tree
  PLLRootedTree &_speciesTree;
  // list of all species tree nodes used for the likelihood computation
  std::vector <corax_rnode_t *> _allSpeciesNodes;
  std::vector <corax_rnode_t *> _prunedSpeciesNodes;
  // map gene leaves to species leaves. The mapping is not computed by this class 
  std::vector<unsigned int> _geneToSpecies;
  // defines at which level (species/genes/none) we do imcremental recomputations  
  PartialLikelihoodMode _likelihoodMode;
  // fraction of missing genes, indexed by gene ids
  std::vector<double> _fm;
  // number of genes that cover each species leaf
  std::vector<unsigned int> _speciesCoverage;
  // number of species leaf that are covered by at least one gene
  unsigned int _numberOfCoveredSpecies;
  std::map<std::string, std::string> _geneNameToSpeciesName;
  // species name (leaf label) to species node index 
  std::map<std::string, unsigned int> _speciesNameToId;
  // if true, we have to recompute the whole likelihood from scratch 
  // (no imcremental recomputation)
  bool _allSpeciesNodesInvalid;
  // species internal nodes that need to be recomputed
  std::unordered_set<corax_rnode_t *> _invalidatedSpeciesNodes;
  //Internal representation of the current species tree. Always use these
  //pointers to be compliant with the pruned species tree mode
  std::vector<corax_rnode_t *> _speciesLeft;
  std::vector<corax_rnode_t *> _speciesRight;
  std::vector<corax_rnode_t *> _speciesParent;
  corax_rnode_t *_prunedRoot;

private:
  size_t getTreeHashRec(const corax_rnode_t *node, size_t i) const;
};

