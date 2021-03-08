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
//#define IS_PROBA(x) true
//#define ASSERT_PROBA(x) assert(true);

double log(ScaledValue v);
typedef std::vector< std::vector <double> > RatesVector;

static bool fillNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd = nullptr)  
{
  bool addMyself = true;
  if (nodesToAdd) {
    addMyself = (nodesToAdd->find(node) != nodesToAdd->end());
  }
  if (node->left) {
    assert(node->right);
    addMyself |= fillNodesPostOrder(node->left, nodes, nodesToAdd);
    addMyself |= fillNodesPostOrder(node->right, nodes, nodesToAdd);
  }
  if (addMyself) {
    nodes.push_back(node);
  }
  return addMyself;
}

static pll_unode_t *getOther(pll_unode_t *ref, pll_unode_t *n1, pll_unode_t *n2)
{
  return (ref == n1) ? n2 : n1;
}

template<class REAL> 
REAL getRandom(REAL max)
{
  return max * Random::getProba();
}
template<class C, class REAL>
int sampleIndex(const C &container)
{
  REAL sum = REAL();
  for (auto value: container) {
    sum += value;
  }
  if (sum == REAL()) {
    return -1;
  }
  REAL stopAt = getRandom(sum);
  sum = REAL();
  unsigned int index = 0;
  for (auto value: container) {
    sum += value;
    if (stopAt < sum) {
      return index;
    }
    index++;
  }
  assert(false);
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
  
  /**
   * (incrementally) compute and return the likelihood of the gene tree 
   */
  virtual double computeLogLikelihood() = 0;


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

  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate) = 0;

  virtual void setFractionMissingGenes(const std::string &fractionMissingFile) = 0;

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
  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);

  pll_rnode_t *getSpeciesLeft(pll_rnode_t *node) {
    return _speciesLeft[node->node_index];
  }

  pll_rnode_t *getSpeciesRight(pll_rnode_t *node) {
    return _speciesRight[node->node_index];
  }

  pll_rnode_t *getSpeciesParent(pll_rnode_t *node) {
    return _speciesParent[node->node_index];
  }

  pll_rnode_t *getPrunedRoot() {return _prunedRoot;}
protected:
  virtual void initSpeciesTree();
  virtual void recomputeSpeciesProbabilities() = 0;
  virtual void beforeComputeLogLikelihood(); 
  
  bool fillPrunedNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd = nullptr);  



protected:
  RecModelInfo _info;
  PLLRootedTree &_speciesTree;
  std::vector <pll_rnode_t *> _speciesNodesToUpdate;
  std::vector <pll_rnode_t *> _allSpeciesNodes;
  unsigned int _allSpeciesNodesCount;
  std::vector<unsigned int> _geneToSpecies;
  PartialLikelihoodMode _likelihoodMode;
  std::vector<double> _fm;
  std::vector<unsigned int> _speciesCoverage;
  unsigned int _numberOfCoveredSpecies;
  std::map<std::string, std::string> _geneNameToSpeciesName;
  std::map<std::string, unsigned int> _speciesNameToId;
  bool _allSpeciesNodesInvalid;
  std::unordered_set<pll_rnode_t *> _invalidatedSpeciesNodes;
  // left, right and parent species vectors, 
  // indexed with the species nodex_index
  std::vector<pll_rnode_t *> _speciesLeft;
  std::vector<pll_rnode_t *> _speciesRight;
  std::vector<pll_rnode_t *> _speciesParent;
  pll_rnode_t *_prunedRoot;

};

