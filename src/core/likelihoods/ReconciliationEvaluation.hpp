#pragma once

#include <enums.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>

#include <memory>

using namespace std;
  
/**
 *  Wrapper around the reconciliation likelihood classes
 */
class ReconciliationEvaluation {
public:
  
  /**
   *  Constructor 
   *  @param speciesTree: rooted species tree (fixed)
   *  @param map: gene-to-species mapping
   *  @param reconciliationModelStr: the reconciliation model to use
   *  @param rootedGeneTree: should we compute the likelihood of a rooted or unrooted gene tree?
   */
  ReconciliationEvaluation(pll_rtree_t *speciesTree,
    const GeneSpeciesMapping& map,
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

  pll_unode_t *getRoot() {return reconciliationModel->getRoot();}
  void setRoot(pll_unode_t * root) {reconciliationModel->setRoot(root);}
  
  bool implementsTransfers() {return reconciliationModel->implementsTransfers();}

  /**
   *  @param input treeinfo
   *  @return the reconciliation likelihood of this tree
   */
  double evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo);

  /**
   *  Invalidate the CLV at a given node index
   *  Must be called on the nodes affected by a move 
   */
  void invalidateCLV(int nodeIndex);
  
  
  void inferMLScenario(shared_ptr<pllmod_treeinfo_t> treeinfo, Scenario &scenario) {
    reconciliationModel->inferMLScenario(treeinfo, scenario);
  }
private:
  shared_ptr<AbstractReconciliationModel> getRecModel(RecModel recModel);
  shared_ptr<AbstractReconciliationModel> reconciliationModel;
};

