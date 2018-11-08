#ifndef _ABSTRACT_MODEL_HPP_
#define _ABSTRACT_MODEL_HPP_

#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>

class AbstractReconciliationModel {
public:
  virtual ~AbstractReconciliationModel() {};
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0) = 0;
  virtual void setSpeciesTree(pll_rtree_t *geneTree) = 0;
  virtual void setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo) = 0;
  virtual void setGeneSpeciesMap(const GeneSpeciesMapping &map) = 0;
  virtual void setRoot(pll_unode_t * root) = 0;
  virtual pll_unode_t *getRoot() = 0;
  virtual double computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo) = 0;
};


#endif
