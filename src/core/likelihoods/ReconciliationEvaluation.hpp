#pragma once

#include <util/enums.hpp>
#include <util/RecModelInfo.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <vector>
#include <memory>
#include <maths/DTLRates.hpp>
#include <maths/Parameters.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <trees/PLLRootedTree.hpp>

class ReconciliationModelInterface;
class Scenario;

/**
 *  Wrapper around the reconciliation likelihood classes
 */
class ReconciliationEvaluation {
public:
  
  /**
   *  Constructor 
   *  @param speciesTree rooted species tree (std::fixed)
   *  @param initialGeneTree initial gene tree
   *  @param geneSpeciesMapping gene-to-species geneSpeciesMapping
   *  @param recModelInfo description of the reconciliation model
   */
  ReconciliationEvaluation(PLLRootedTree &speciesTree,
    PLLUnrootedTree &initialGeneTree,
    const GeneSpeciesMapping& geneSpeciesMapping,
    const RecModelInfo &recModelInfo); 
  
  /**
   * Forbid copy
   */
  ReconciliationEvaluation(const ReconciliationEvaluation &) = delete;
  ReconciliationEvaluation & operator = (const ReconciliationEvaluation &) = delete;
  ReconciliationEvaluation(ReconciliationEvaluation &&) = delete;
  ReconciliationEvaluation & operator = (ReconciliationEvaluation &&) = delete;
  ~ReconciliationEvaluation();
  void setRates(const Parameters &parameters);

  /**
   * Get the current root of the gene tree. Return null if the tree does not have a 
   * current root (in unrooted mode)
   * This method is mostly used for rollbacking to a previous state. In most of the
   * cases you should call inferMLRoot instead.
   */
  pll_unode_t *getRoot();
  void setRoot(pll_unode_t * root);

  void enableMADRooting(bool enable);

  /**
   *  @param input geneTree
   */
  double evaluate();

  bool implementsTransfers() {return Enums::accountsForTransfers(_recModelInfo.model);} 

  /*
   *  Call this everytime that the species tree changes
   */
  void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);

  void setPartialLikelihoodMode(PartialLikelihoodMode mode);

  /**
   *  Invalidate the CLV at a given node index
   *  Must be called on the nodes affected by a move 
   */
  void invalidateCLV(unsigned int nodeIndex);
 
  void invalidateAllCLVs();
  void invalidateAllSpeciesCLVs();
 
  pll_unode_t *inferMLRoot();
  
  void inferMLScenario(Scenario &scenario, bool stochastic = false);

  RecModel getRecModel() const {return _recModelInfo.model;}
  
private:
  PLLRootedTree &_speciesTree;
  PLLUnrootedTree &_initialGeneTree;
  GeneSpeciesMapping _geneSpeciesMapping;
  RecModelInfo _recModelInfo; 
  bool _infinitePrecision;
  std::vector<std::vector<double> > _rates;
  // we actually own this pointer, but we do not 
  // wrap it into a unique_ptr to allow forward definition
  ReconciliationModelInterface *_evaluators;
private:
  ReconciliationModelInterface *buildRecModelObject(RecModel recModel, bool infinitePrecision);
  pll_unode_t *computeMLRoot();
  void updatePrecision(bool infinitePrecision);
};
  
using Evaluations = std::vector<std::shared_ptr<ReconciliationEvaluation> >;
using PerCoreEvaluations = std::vector<std::shared_ptr<ReconciliationEvaluation> >;

