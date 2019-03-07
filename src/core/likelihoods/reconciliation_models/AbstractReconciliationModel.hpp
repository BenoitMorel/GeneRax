#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
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
  virtual void init(pll_rtree_t *speciesTree, const GeneSpeciesMapping &map, bool rootedGeneTree);
  
  virtual ~AbstractReconciliationModel() {};
 
  /**
   * Set the DTL rates, and update probabilities relative to the species tree only
   */
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0) = 0;
  
  /**
   *  Does this model support transfer? Relevant to know which parameter to optimize
   */
  virtual bool implementsTransfers() = 0;
 
  /**
   * (incrementally) compute and return the likelihood of the input gene tree 
   */
  double computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  
  
  /**
   *  Get/set the root of the gene tree (only relevant in rooted gene tree mode)
   */ 
  void setRoot(pll_unode_t * root) {geneRoot_ = root;}
  pll_unode_t *getRoot() {return geneRoot_;}
 
  /**
   * invalidate one or all the gene CLVs
   */
  void invalidateAllCLVs();
  void invalidateCLV(int geneNodeIndex);

  /**
   *  Fill scenario with the maximum likelihood set of 
   *  events that would lead to treeinfo
   **/
  void inferMLScenario(shared_ptr<pllmod_treeinfo_t> treeinfo, Scenario &scenario);

protected:
  // called by the constructor
  virtual void setSpeciesTree(pll_rtree_t *speciesTree);

  // Called when computeLogLikelihood is called for the first time
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  // Called by computeLogLikelihood
  virtual void updateCLV(pll_unode_t *geneNode) = 0;
  // Called by computeLogLikelihood
  virtual void computeRootLikelihood(pllmod_treeinfo_t &treeinfo,
    pll_unode_t *virtualRoot) = 0;
  // Called by computeLogLikelihood
  virtual ScaledValue getRootLikelihood(pllmod_treeinfo_t &treeinfo,
    pll_unode_t *root) const = 0;
  virtual ScaledValue getRootLikelihood(pllmod_treeinfo_t &treeinfo,
    pll_unode_t *root, pll_rnode_t *speciesRoot) = 0;
  // Called by inferMLScenario
  // fills scenario with the best likelihood set of events that 
  // would lead to the subtree of geneNode under speciesNode
  // Can assume that all the CLVs are filled
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false) = 0;
 
  
  void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds);
  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  void updateRoot(pllmod_treeinfo_t &treeinfo);
  
  /**
   *  - In rooted gene tree mode, and if the gene tree already has a virtual root,
   *  return this root and its direct neighbors
   *  - Else, return all the possible virtual roots
   */
  void getRoots(pllmod_treeinfo_t &treeinfo, 
    vector<pll_unode_t *> &roots,
    const vector<int> &geneIds);

  /**
   *  Get the left or right child of node. If node is a virtual 
   *  root, the implementation is different
   */
  static pll_unode_t *getLeft(pll_unode_t *node, bool virtualRoot);  
  static pll_unode_t *getRight(pll_unode_t *node, bool virtualRoot) ;  
  
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
protected:
  bool rootedGeneTree_;
  pll_unode_t *geneRoot_;
  int speciesNodesCount_;
  vector <pll_rnode_t *> speciesNodes_;
  pll_rtree_t *speciesTree_;
  vector<int> geneToSpecies_;
  bool firstCall_;
  // gene ids in postorder 
  vector<int> _geneIds;
  int _maxGeneId;

private:
  pll_unode_t *computeMLRoot(pllmod_treeinfo_t &treeinfo);
  void computeMLRoot(pllmod_treeinfo_t &treeinfo, 
    pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot);
  virtual void computeLikelihoods(pllmod_treeinfo_t &treeinfo);
  double getSumLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  void updateCLVsRec(pll_unode_t *node);
  void markInvalidatedNodes(pllmod_treeinfo_t &treeinfo);
  void markInvalidatedNodesRec(pll_unode_t *node);
  void fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) ;
  map<string, string> geneNameToSpeciesName_;
  map<string, int> speciesNameToId_;
  
  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  unordered_set<int> _invalidatedNodes;

  // is the CLV up to date?
  vector<bool> _isCLVUpdated;

};

