#pragma once

#include <enums.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <vector>
#include <memory>
#include <maths/DTLRates.hpp>

  
/**
 *  Wrapper around the reconciliation likelihood classes
 */
class ReconciliationEvaluation {
public:
  
  /**
   *  Constructor 
   *  @param speciesTree rooted species tree (std::fixed)
   *  @param geneSpeciesMapping gene-to-species geneSpeciesMappingping
   *  @param reconciliationModelStr the reconciliation model to use
   *  @param rootedGeneTree should we compute the likelihood of a rooted or unrooted gene tree?
   */
  ReconciliationEvaluation(pll_rtree_t *speciesTree,
    const GeneSpeciesMapping& geneSpeciesMapping,
    RecModel recModel,
    bool rootedGeneTree);

  /**
   *  Set the DTL rates
   *  @param dupRate
   *  @param lossRate
   *  @param transferRate
   */ 
  void setRates(double dupRate, double lossRate, 
    double transferRate = 0.0);

  /*
   * Set the per-species lineage rates
   *  @param ratesVector DTL rates indexed with species node index
   */
  void setRates(const DTLRatesVector &ratesVector);

  AbstractReconciliationModel *getReconciliationModel() {return reconciliationModel.get();}

  pll_unode_t *getRoot() {return reconciliationModel->getRoot();}
  void setRoot(pll_unode_t * root) {reconciliationModel->setRoot(root);}

  double evaluate(pll_utree_t *utree);
  /**
   *  @param input treeinfo
   *  @return the reconciliation likelihood of this tree
   */
  double evaluate(std::shared_ptr<pllmod_treeinfo_t> treeinfo);

  bool implementsTransfers() {return Enums::accountsForTransfers(_model);} 

  /**
   *  Invalidate the CLV at a given node index
   *  Must be called on the nodes affected by a move 
   */
  void invalidateCLV(unsigned int nodeIndex);
  
  
  void inferMLScenario(Scenario &scenario) {
    reconciliationModel->inferMLScenario(scenario);
  }

private:
  std::shared_ptr<AbstractReconciliationModel> getRecModelObject(RecModel recModel);
  std::shared_ptr<AbstractReconciliationModel> reconciliationModel;
  RecModel _model;
  unsigned int _speciesCount;
};

