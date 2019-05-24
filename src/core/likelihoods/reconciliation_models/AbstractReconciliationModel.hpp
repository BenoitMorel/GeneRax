#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/SubtreeRepeatsCache.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <Scenario.hpp>

#include <unordered_set>
#include <maths/ScaledValue.hpp>

/**
 *  Interface and common implementations for 
 *  all the reconciliation likelihood computation
 *  classes
 */
class AbstractReconciliationModel {
public:
  AbstractReconciliationModel();
     
  /**
   *  Has to be called just after the constructor before anything else
   *  We do not call it in the constructor because it is virtual and
   *  calls virtual functions
   */
  virtual void init(pll_rtree_t *speciesTree, const GeneSpeciesMapping &geneSpeciesMappingp, bool rootedGeneTree);
  
  virtual ~AbstractReconciliationModel() {};
 
  /**
   * Set the DTL rates, and update probabilities relative to the species tree only
   */
  void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  
  /*
   * Set the per-species lineage rates
   *  @param dupRate
   *  @param lossRate
   *  @param transferRate
   */
  virtual void setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates) = 0;

  
  /**
   * (incrementally) compute and return the likelihood of the input gene tree 
   */
  double computeLogLikelihood(pll_utree_t *tree);
  
  
  /**
   *  Get/set the root of the gene tree (only relevant in rooted gene tree mode)
   */ 
  void setRoot(pll_unode_t * root) {geneRoot_ = root;}
  pll_unode_t *getRoot() {return geneRoot_;}
 
  /**
   * invalidate one or all the gene CLVs
   */
  void invalidateAllCLVs();
  void invalidateCLV(unsigned int geneNodeIndex);

  /**
   *  Fill scenario with the maximum likelihood set of 
   *  events that would lead to the  current tree
   **/
  void inferMLScenario(Scenario &scenario);
  
  
  SubtreeRepeatsCache &getCache() {return _cache;}

protected:
  // called by the constructor
  virtual void setSpeciesTree(pll_rtree_t *speciesTree);

  // Called when computeLogLikelihood is called for the first time
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // Called by computeLogLikelihood
  virtual void updateCLV(pll_unode_t *geneNode) = 0;
  // Called by computeLogLikelihood
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot) = 0;
  // Called by computeLogLikelihood
  virtual ScaledValue getRootLikelihood(pll_unode_t *root) const = 0;
  virtual ScaledValue getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) = 0;
  // Called by inferMLScenario
  // fills scenario with the best likelihood set of events that 
  // would lead to the subtree of geneNode under speciesNode
  // Can assume that all the CLVs are filled
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false) = 0;
 
  
  void initFromUtree(pll_utree_t *tree);
  /**
   *  - In rooted gene tree mode, and if the gene tree already has a virtual root,
   *  return this root and its direct neighbors
   *  - Else, return all the possible virtual roots
   */
  void getRoots(std::vector<pll_unode_t *> &roots,
    const std::vector<unsigned int> &geneIds);

  /**
   *  Get the left or right child of node. If node is a virtual 
   *  root, the implementation is different
   */
  pll_unode_t *getLeft(pll_unode_t *node, bool virtualRoot) const;  
  pll_unode_t *getRight(pll_unode_t *node, bool virtualRoot) const; 
  pll_unode_t *getLeftRepeats(pll_unode_t *node, bool virtualRoot);  
  pll_unode_t *getRightRepeats(pll_unode_t *node, bool virtualRoot) ; 

  
  void updateCLVs();
protected:
  bool rootedGeneTree_;
  pll_unode_t *geneRoot_;
  unsigned int speciesNodesCount_;
  std::vector <pll_rnode_t *> speciesNodes_;
  pll_rtree_t *speciesTree_;
  std::vector<unsigned int> geneToSpecies_;
  bool firstCall_;
  // gene ids in postorder 
  std::vector<unsigned int> _geneIds;
  unsigned int _maxGeneId;

private:
  void mapGenesToSpecies();
  pll_unode_t *computeMLRoot();
  void computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot);
  virtual void computeLikelihoods();
  double getSumLikelihood();
  void updateCLVsRec(pll_unode_t *node);
  void markInvalidatedNodes();
  void markInvalidatedNodesRec(pll_unode_t *node);
  void fillNodesPostOrder(pll_rnode_t *node, std::vector<pll_rnode_t *> &nodes) ;
  std::map<std::string, std::string> geneNameToSpeciesName_;
  std::map<std::string, unsigned int> speciesNameToId_;
  
  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  std::unordered_set<unsigned int> _invalidatedNodes;

  // is the CLV up to date?
  std::vector<bool> _isCLVUpdated;
  std::vector<pll_unode_t *> _allNodes;
  SubtreeRepeatsCache _cache;
};


