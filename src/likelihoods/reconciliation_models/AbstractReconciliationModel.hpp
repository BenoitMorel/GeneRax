#ifndef _ABSTRACT_MODEL_HPP_
#define _ABSTRACT_MODEL_HPP_

#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>

class AbstractReconciliationModel {
public:
  AbstractReconciliationModel();
  virtual ~AbstractReconciliationModel() {};
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0) = 0;
  virtual double computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo) = 0;
  
  virtual void setSpeciesTree(pll_rtree_t *geneTree);
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo);
  virtual void setGeneSpeciesMap(const GeneSpeciesMapping &map);
  virtual void setRoot(pll_unode_t * root) {geneRoot_ = root;}
  virtual pll_unode_t *getRoot() {return geneRoot_;}

protected:
  void getIdsPostOrder(pllmod_treeinfo_t &tree, vector<int> &nodeIds);
  void mapGenesToSpecies(pllmod_treeinfo_t &treeinfo);
  
protected:
  pll_unode_t *geneRoot_;
  int speciesNodesCount_;
  vector <pll_rnode_t *> speciesNodes_;
  pll_rtree_t *speciesTree_;
  vector<int> geneToSpecies_;

private:
  void fillNodesPostOrder(pll_rnode_t *node, vector<pll_rnode_t *> &nodes) ;
  map<string, string> geneNameToSpeciesName_;
  map<string, int> speciesNameToId_;

};


#endif
