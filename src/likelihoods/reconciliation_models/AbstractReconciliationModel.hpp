#ifndef _ABSTRACT_MODEL_HPP_
#define _ABSTRACT_MODEL_HPP_

#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>

class AbstractReconciliationModel {
public:
  AbstractReconciliationModel();
  virtual ~AbstractReconciliationModel() {};
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0) = 0;
  virtual double computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo) = 0;
  
  virtual void setSpeciesTree(pll_rtree_t *geneTree);
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual void setGeneSpeciesMap(const GeneSpeciesMapping &map);
  virtual void setRoot(pll_unode_t * root) {geneRoot_ = root;}
  virtual void invalidateCLV(int nodeIndex) {}
  virtual pll_unode_t *getRoot() {return geneRoot_;}

protected:
  void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds);
  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  void getRoots(pllmod_treeinfo_t &treeinfo, 
    vector<pll_unode_t *> &roots,
    const vector<int> &geneIds);
  size_t getPrunedSpeciesCount() {return prunedSpeciesNodes_.size();}
  int getPrunedSpeciesIndex(pll_rnode_t *node) {return speciesToPrunedSpecies_[node->node_index];}
  const vector<pll_rnode_t *> &getPrunedSpecies() {return prunedSpeciesNodes_;}
  pll_unode_t *getLeft(pll_unode_t *node, bool virtualRoot);  
  pll_unode_t *getRight(pll_unode_t *node, bool virtualRoot);  
protected:
  pll_unode_t *geneRoot_;
  int speciesNodesCount_;
  vector<pll_rnode_t *> speciesNodes_;
  pll_rtree_t *speciesTree_;
  vector<int> geneToSpecies_;
  
  vector<int> speciesToPrunedSpecies_;
  vector<pll_rnode_t *> prunedSpeciesNodes_;

private:
  void fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) ;
  map<string, string> geneNameToSpeciesName_;
  map<string, int> speciesNameToId_;

};


#endif
