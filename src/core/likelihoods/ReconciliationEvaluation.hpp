#pragma once

#include <util/enums.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <vector>
#include <memory>
#include <maths/DTLRates.hpp>
#include <maths/Parameters.hpp>

double log(ScaledValue v) ;
  
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

  void setRates(const Parameters &parameters);

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


  pll_unode_t *getRoot() {return _reconciliationModel->getRoot();}
  void setRoot(pll_unode_t * root) {_reconciliationModel->setRoot(root);}


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
 
  pll_unode_t *inferMLRoot(pll_utree_t *utree);
  
  void inferMLScenario(pll_utree_t *tree, Scenario &scenario);

  RecModel getRecModel() const {return _model;}
private:
  pll_unode_t *computeMLRoot();
  pll_rtree_t *_speciesTree;
  unsigned int _speciesCount;
  GeneSpeciesMapping _geneSpeciesMapping;
  bool _rootedGeneTree;
  RecModel _model; 
  bool _infinitePrecision;
  std::vector<double> _dupRates;
  std::vector<double> _lossRates;
  std::vector<double> _transferRates;
  std::unique_ptr<ReconciliationModelInterface> _reconciliationModel;
  std::unique_ptr<ReconciliationModelInterface> buildRecModelObject(RecModel recModel, bool infinitePrecision);
  void updatePrecision(bool infinitePrecision);
};

