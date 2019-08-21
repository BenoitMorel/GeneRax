#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>


/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In addition, we forbid transfers to parent species
*/
class UndatedDTLModel: public AbstractReconciliationModel {
public:
  UndatedDTLModel();
  virtual ~UndatedDTLModel();
  
  // overloaded from parent
  virtual void setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates);
protected:
  // overloaded from parent
  virtual void setInitialGeneTree(pll_utree_t *tree);
  // overloaded from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overloaded from parent
  virtual double getRootLikelihood(pll_unode_t *root) const;
  // overload from parent
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot);
  virtual double getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {return _uq[root->node_index + _maxGeneId + 1][speciesRoot->node_index];}
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false);
private:
  // model
  std::vector<double> _PD; // Duplication probability, per branch
  std::vector<double> _PL; // Loss probability, per branch
  std::vector<double> _PT; // Transfer probability, per branch
  std::vector<double> _PS; // Speciation probability, per branch

  // SPECIES
  std::vector<double> _uE; // Probability for a gene to become extinct on each brance
  double _transferExtinctionSum;
  std::vector<double> _ancestralExctinctionCorrection;  

  // CLVs
  // _uq[geneId][speciesId] = probability of a gene node rooted at a species node
  // to produce the subtree of this gene node
  std::vector<std::vector<double> > _uq;
  std::vector<double> _survivingTransferSums;
  std::vector<std::vector<double> > _ancestralCorrection;



private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      double &proba,
      bool isVirtualRoot = false);
  void updateTransferSums(double &transferExtinctionSum,
    std::vector<double> &ancestralExtinctionCorrection,
    const std::vector<double> &probabilities);
  void resetTransferSums(double &transferSum,
    std::vector<double> &ancestralCorrection);
  
  double getCorrectedTransferExtinctionSum(unsigned int speciesNode) const {
  return _transferExtinctionSum - _ancestralExctinctionCorrection[speciesNode];
  }

  double getCorrectedTransferSum(unsigned int geneId, unsigned int speciesId) const
  {
    return _survivingTransferSums[geneId] - _ancestralCorrection[geneId][speciesId];

  }

  pll_rnode_t *getBestTransfer(unsigned int gid, pll_rnode_t *speciesNode); 
};

