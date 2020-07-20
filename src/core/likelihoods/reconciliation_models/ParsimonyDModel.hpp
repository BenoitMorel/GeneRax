#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <util/Scenario.hpp>
#include <cmath>





class ParsimonyDModel: public AbstractReconciliationModel<double> {
public:
  ParsimonyDModel(PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMappingp, 
      bool rootedGeneTree,
      double minGeneBranchLength,
      bool pruneSpeciesTree):
    AbstractReconciliationModel<double>(speciesTree, geneSpeciesMappingp, rootedGeneTree, minGeneBranchLength, pruneSpeciesTree),
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
  virtual void updateCLV(pll_unode_t *geneNode);
  // overload from parent
  virtual double getGeneRootLikelihood(pll_unode_t *root) const;
  virtual double  getGeneRootLikelihood(pll_unode_t *root, pll_rnode_t *) {
    return _dlclvs[root->node_index + this->_maxGeneId + 1].cost;
  }

  // overload from parent
  virtual void recomputeSpeciesProbabilities() {}
  virtual double getLikelihoodFactor() const {return 1.0;}
  // overload from parent
  virtual void computeGeneRootLikelihood(pll_unode_t *virtualRoot);
  // overlead from parent
  virtual void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      double &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr,
      bool stochastic = false);
  virtual bool isParsimony() const {return true;}
  virtual bool sumOverAllOriginations() const {return false;}
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
  std::vector<pll_rnode_s *> &getSpeciesNodesToUpdate() {
    return this->_speciesNodesToUpdate;
  }

};



void ParsimonyDModel::setInitialGeneTree(PLLUnrootedTree &tree)
{
  AbstractReconciliationModel<double>::setInitialGeneTree(tree);
  _dlclvs = std::vector<DLCLV>(2 * (this->_maxGeneId + 1));
}


void ParsimonyDModel::updateCLV(pll_unode_t *geneNode)
{
  computeProbability(geneNode, 
        nullptr, 
        _dlclvs[geneNode->node_index].cost);
}

static void getPolytomyD(pll_unode_t *node, 
    std::vector<pll_unode_t *> &polytomy,
    double minBranchLength = 0.0000011)
{
  if (!node->next || node->length > minBranchLength) {
    polytomy.push_back(node);
  } else {
    getPolytomyD(node->next->back, polytomy);
    getPolytomyD(node->next->next->back, polytomy);
  }
}


void ParsimonyDModel::computeProbability(pll_unode_t *geneNode, 
    pll_rnode_t *, 
      double &proba,
      bool isVirtualRoot,
      Scenario *,
      Scenario::Event *,
      bool)
  
{
  proba = 0.0;
  auto gid = geneNode->node_index;
  if (!geneNode->next) {
    return; 
  }
  auto leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
  auto rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  auto v = leftGeneNode->node_index;
  auto w = rightGeneNode->node_index;
  proba = _dlclvs[v].cost + _dlclvs[w].cost;
  std::vector<pll_unode_t *> polytomy;
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
  

double ParsimonyDModel::getGeneRootLikelihood(pll_unode_t *root) const
{
  auto u = root->node_index + this->_maxGeneId + 1;
  return _dlclvs[u].cost;
}


void ParsimonyDModel::computeGeneRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  computeProbability(virtualRoot, nullptr, _dlclvs[u].cost, true);
}




