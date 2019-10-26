#pragma once

#include <util/enums.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <vector>
#include <memory>
#include <maths/DTLRates.hpp>
#include <maths/Parameters.hpp>
#include <trees/PLLUnrootedTree.hpp>

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
  
  /**
   * Forbid copy
   */
  ReconciliationEvaluation(const ReconciliationEvaluation &) = delete;
  ReconciliationEvaluation & operator = (const ReconciliationEvaluation &) = delete;
  ReconciliationEvaluation(ReconciliationEvaluation &&) = delete;
  ReconciliationEvaluation & operator = (ReconciliationEvaluation &&) = delete;

  void setTransferFrequencies(const Parameters &parameters);
  void setRates(const Parameters &parameters);


  /**
   * Get the current root of the gene tree. Return null if the tree does not have a 
   * current root (in unrooted mode)
   * This method is mostly used for rollbacking to a previous state. In most of the
   * cases you should call inferMLRoot instead.
   */
  pll_unode_t *getRoot() {return _reconciliationModel->getRoot();}
  void setRoot(pll_unode_t * root) {_reconciliationModel->setRoot(root);}

  /**
   *  @param input utree
   *  @return the reconciliation likelihood of this tree
   */
  double evaluate(PLLUnrootedTree &geneTree);

  bool implementsTransfers() {return Enums::accountsForTransfers(_model);} 

  /**
   *  Invalidate the CLV at a given node index
   *  Must be called on the nodes affected by a move 
   */
  void invalidateCLV(unsigned int nodeIndex);
 
  pll_unode_t *inferMLRoot(PLLUnrootedTree &geneTree);
  
  void inferMLScenario(PLLUnrootedTree &tree, Scenario &scenario);

  RecModel getRecModel() const {return _model;}
private:
  pll_rtree_t *_speciesTree; // I do not own this one
  unsigned int _speciesCount;
  GeneSpeciesMapping _geneSpeciesMapping;
  bool _rootedGeneTree;
  RecModel _model; 
  bool _infinitePrecision;
  std::vector<double> _dupRates;
  std::vector<double> _lossRates;
  std::vector<double> _transferRates;
  std::vector< std::vector<double> > _transferFrequencies;
  std::unique_ptr<ReconciliationModelInterface> _reconciliationModel;
private:
  
  std::unique_ptr<ReconciliationModelInterface> buildRecModelObject(RecModel recModel, bool infinitePrecision);
  pll_unode_t *computeMLRoot();
  void updatePrecision(bool infinitePrecision);
};

