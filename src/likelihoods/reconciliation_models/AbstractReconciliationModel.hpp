#ifndef _ABSTRACT_MODEL_HPP_
#define _ABSTRACT_MODEL_HPP_

#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>


/**
 *  Interface and common implementations for 
 *  all the reconciliation likelihood computation
 *  classes
 */
class AbstractReconciliationModel {
public:
  AbstractReconciliationModel(pll_rtree_t *speciesTree, const GeneSpeciesMapping &map);
  
  virtual ~AbstractReconciliationModel() {};
  
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0) = 0;
  
  virtual double computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  
  virtual void invalidateCLV(int nodeIndex) {}
  

protected:
  // called by the constructor
  virtual void setSpeciesTree(pll_rtree_t *geneTree);
  virtual void setGeneSpeciesMap(const GeneSpeciesMapping &map);

  // Called when computeLogLikelihood is called for the first time
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  // Called by computeLogLikelihood
  virtual double computeLogLikelihoodInternal(shared_ptr<pllmod_treeinfo_t> treeinfo) = 0;
 
  virtual void setRoot(pll_unode_t * root) {geneRoot_ = root;}

  virtual pll_unode_t *getRoot() {return geneRoot_;}
  
  void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds);
  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  void getRoots(pllmod_treeinfo_t &treeinfo, 
    vector<pll_unode_t *> &roots,
    const vector<int> &geneIds);
  pll_unode_t *getLeft(pll_unode_t *node, bool virtualRoot);  
  pll_unode_t *getRight(pll_unode_t *node, bool virtualRoot);  
protected:
  pll_unode_t *geneRoot_;
  int speciesNodesCount_;
  vector <pll_rnode_t *> speciesNodes_;
  pll_rtree_t *speciesTree_;
  vector<int> geneToSpecies_;
  bool firstCall_;

private:
  void fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) ;
  map<string, string> geneNameToSpeciesName_;
  map<string, int> speciesNameToId_;

};


#endif
