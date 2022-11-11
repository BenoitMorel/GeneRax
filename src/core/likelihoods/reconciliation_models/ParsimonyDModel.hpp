#pragma once

#include <likelihoods/reconciliation_models/GTBaseReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <util/Scenario.hpp>
#include <cmath>





class ParsimonyDModel: public GTBaseReconciliationModel<double> {
public:
  ParsimonyDModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMappingp, 
      const RecModelInfo &recModelInfo):
    GTBaseReconciliationModel<double>(speciesTree, 
        geneSpeciesMappingp, 
        recModelInfo),
    _costD(-1.0)
    {}
  
  
  ParsimonyDModel(const ParsimonyDModel &) = delete;
  ParsimonyDModel & operator = (const ParsimonyDModel &) = delete;
  ParsimonyDModel(ParsimonyDModel &&) = delete;
  ParsimonyDModel & operator = (ParsimonyDModel &&) = delete;
  virtual ~ParsimonyDModel() {}
  
  // overloaded from parent
  virtual void setRates(const RatesVector &){}
protected:
  // overload from parent
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  // overload from parent
  virtual void updateCLV(corax_unode_t *geneNode);
  // overload from parent
  virtual double getGeneRootLikelihood(corax_unode_t *root) const;
  virtual double  getGeneRootLikelihood(corax_unode_t *root, corax_rnode_t *) {
    return _dlclvs[root->node_index + this->_maxGeneId + 1].cost;
  }
  // overload from parent
  virtual void recomputeSpeciesProbabilities() {}
  virtual double getLikelihoodFactor() const {return 1.0;}
  // overload from parent
  virtual void computeGeneRootLikelihood(corax_unode_t *virtualRoot);
  // overlead from parent
  virtual void computeProbability(corax_unode_t *geneNode, corax_rnode_t *speciesNode, 
      double &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr,
      bool stochastic = false);
  virtual bool isParsimony() const {return true;}
private:
  // parsimony costs
  double _costD;
  
  // uq[geneId][speciesId] = parsimony cost
  // of a gene node rooted at a species node
  // to produce the subtree of this gene node
  struct DLCLV {
    double cost;
    DLCLV():  cost(-std::numeric_limits<double>::infinity()) {}
  };
  std::vector<DLCLV> _dlclvs;
 
private:

};



void ParsimonyDModel::setInitialGeneTree(PLLUnrootedTree &tree)
{
  GTBaseReconciliationModel<double>::setInitialGeneTree(tree);
  _dlclvs = std::vector<DLCLV>(2 * (this->_maxGeneId + 1));
}


void ParsimonyDModel::updateCLV(corax_unode_t *geneNode)
{
  computeProbability(geneNode, 
        nullptr, 
        _dlclvs[geneNode->node_index].cost);
}

static void getPolytomyD(corax_unode_t *node, 
    std::vector<corax_unode_t *> &polytomy,
    double minBranchLength = 0.0000011)
{
  if (!node->next || node->length > minBranchLength) {
    polytomy.push_back(node);
  } else {
    getPolytomyD(node->next->back, polytomy);
    getPolytomyD(node->next->next->back, polytomy);
  }
}


void ParsimonyDModel::computeProbability(corax_unode_t *geneNode, 
    corax_rnode_t *, 
      double &proba,
      bool isVirtualRoot,
      Scenario *,
      Scenario::Event *,
      bool)
  
{
  proba = 0.0;
  if (!geneNode->next) {
    return; 
  }
  auto leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
  auto rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  auto v = leftGeneNode->node_index;
  auto w = rightGeneNode->node_index;
  proba = _dlclvs[v].cost + _dlclvs[w].cost;
  std::vector<corax_unode_t *> polytomy;
  getPolytomyD(leftGeneNode, polytomy);
  getPolytomyD(rightGeneNode, polytomy);
  if (polytomy.size() > 2) {
    proba = 0.0;
    for (auto node: polytomy) {
      proba += _dlclvs[node->node_index].cost;
    }
    return;
  }
  // D event
  auto leftLCA = this->_geneToSpeciesLCA[v];
  auto rightLCA = this->_geneToSpeciesLCA[w];
  auto &parentsCache = this->_speciesTree.getParentsCache(rightLCA);
  if (parentsCache[leftLCA->node_index]) {
    proba += _costD;
  }
  // else: SL or S event
}
  

double ParsimonyDModel::getGeneRootLikelihood(corax_unode_t *root) const
{
  auto u = root->node_index + this->_maxGeneId + 1;
  return _dlclvs[u].cost;
}


void ParsimonyDModel::computeGeneRootLikelihood(corax_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  computeProbability(virtualRoot, nullptr, _dlclvs[u].cost, true);
}




