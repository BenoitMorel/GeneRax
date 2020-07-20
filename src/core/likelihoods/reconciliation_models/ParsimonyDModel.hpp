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
  virtual double  getGeneRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) {
    return _dlclvs[root->node_index + this->_maxGeneId + 1][speciesRoot->node_index];
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
private:
  // parsimony costs
  double _costD;
  
  // uq[geneId][speciesId] = parsimony cost
  // of a gene node rooted at a species node
  // to produce the subtree of this gene node
  typedef std::vector<double> DLCLV;
  std::vector<DLCLV> _dlclvs;
 
private:
  std::vector<pll_rnode_s *> &getSpeciesNodesToUpdate() {
    return this->_speciesNodesToUpdate;
  }

};



void ParsimonyDModel::setInitialGeneTree(PLLUnrootedTree &tree)
{
  AbstractReconciliationModel<double>::setInitialGeneTree(tree);
  assert(this->_allSpeciesNodesCount);
  assert(this->_maxGeneId);
  std::vector<double> zeros(this->_allSpeciesNodesCount);
  _dlclvs = std::vector<std::vector<double> >(2 * (this->_maxGeneId + 1),zeros);
}


void ParsimonyDModel::updateCLV(pll_unode_t *geneNode)
{
  assert(geneNode);
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    computeProbability(geneNode, 
        speciesNode, 
        _dlclvs[geneNode->node_index][speciesNode->node_index]);
  }
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


void ParsimonyDModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      double &proba,
      bool isVirtualRoot,
      Scenario *,
      Scenario::Event *event,
      bool)
  
{
  assert(!event);
  auto gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  auto e = speciesNode->node_index;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
  }


  proba = -std::numeric_limits<double>::infinity();
  if (isSpeciesLeaf and isGeneLeaf) {
    // present
    if (e == this->_geneToSpecies[gid]) {
      proba = 0.0;
    }
    return;
  }
  if (not isGeneLeaf) {
    auto u_left = leftGeneNode->node_index;
    auto u_right = rightGeneNode->node_index;

    std::vector<pll_unode_t *> polytomy;
    getPolytomyD(leftGeneNode, polytomy);
    getPolytomyD(rightGeneNode, polytomy);
    if (polytomy.size() > 2) {
      proba = 0.0;
      for (auto node: polytomy) {
        proba += _dlclvs[node->node_index][e];
      }
      return;
    }
    // D event
    proba = _costD + _dlclvs[u_left][e] + _dlclvs[u_right][e];
  }
}
  

double ParsimonyDModel::getGeneRootLikelihood(pll_unode_t *root) const
{
  double max =  -std::numeric_limits<double>::infinity();
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    if (max < _dlclvs[u][e]) {
      max = _dlclvs[u][e];
    }
  }
  return max;
}


void ParsimonyDModel::computeGeneRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    auto e = speciesNode->node_index;
    computeProbability(virtualRoot, speciesNode, _dlclvs[u][e], true);
  }
}




