#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <util/Scenario.hpp>
#include <cmath>





class ParsimonyDLModel: public AbstractReconciliationModel<double> {
public:
  ParsimonyDLModel(PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMappingp, 
      bool rootedGeneTree,
      double minGeneBranchLength,
      bool pruneSpeciesTree):
    AbstractReconciliationModel<double>(speciesTree, geneSpeciesMappingp, rootedGeneTree, minGeneBranchLength, pruneSpeciesTree),
    _costS(-0.0),
    _costD(-1.0),
    _costL(-0.0)
    {}
  
  
  ParsimonyDLModel(const ParsimonyDLModel &) = delete;
  ParsimonyDLModel & operator = (const ParsimonyDLModel &) = delete;
  ParsimonyDLModel(ParsimonyDLModel &&) = delete;
  ParsimonyDLModel & operator = (ParsimonyDLModel &&) = delete;
  virtual ~ParsimonyDLModel() {}
  
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
  double _costS;
  double _costD;
  double _costL;
  
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



void ParsimonyDLModel::setInitialGeneTree(PLLUnrootedTree &tree)
{
  AbstractReconciliationModel<double>::setInitialGeneTree(tree);
  assert(this->_allSpeciesNodesCount);
  assert(this->_maxGeneId);
  std::vector<double> zeros(this->_allSpeciesNodesCount);
  _dlclvs = std::vector<std::vector<double> >(2 * (this->_maxGeneId + 1),zeros);
}


void ParsimonyDLModel::updateCLV(pll_unode_t *geneNode)
{
  assert(geneNode);
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    computeProbability(geneNode, 
        speciesNode, 
        _dlclvs[geneNode->node_index][speciesNode->node_index]);
  }
}

static void getPolytomy(pll_unode_t *node, 
    std::vector<pll_unode_t *> &polytomy,
    double minBranchLength = 0.0000011)
{
  if (!node->next || node->length > minBranchLength) {
    polytomy.push_back(node);
  } else {
    getPolytomy(node->next->back, polytomy);
    getPolytomy(node->next->next->back, polytomy);
  }
}


void ParsimonyDLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      double &proba,
      bool isVirtualRoot,
      Scenario *,
      Scenario::Event *event,
      bool)
  
{
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

  if (event) {
    event->geneNode = gid; 
    event->speciesNode = e;
    event->type = ReconciliationEventType::EVENT_None; 
  }


  proba = -std::numeric_limits<double>::infinity();
  if (isSpeciesLeaf and isGeneLeaf) {
    // present
    if (e == this->_geneToSpecies[gid]) {
      proba = 0.0;
    }
    return;
  }
  typedef std::array<double, 5>  ValuesArray;
  ValuesArray values;
  for (unsigned int i = 0; i < 5; ++i) {
    values[i] = -std::numeric_limits<double>::infinity();
  }
 

  if (not isGeneLeaf) {
    auto u_left = leftGeneNode->node_index;
    auto u_right = rightGeneNode->node_index;

    std::vector<pll_unode_t *> polytomy;
    getPolytomy(leftGeneNode, polytomy);
    getPolytomy(rightGeneNode, polytomy);
    if (polytomy.size() > 2) {
      proba = 0.0;
      for (auto node: polytomy) {
        proba += _dlclvs[node->node_index][e];
      }
      return;
    }
    if (not isSpeciesLeaf) {
      // S event
      values[0] = _costS + _dlclvs[u_left][f] + _dlclvs[u_right][g];
      values[1] = _costS + _dlclvs[u_left][g] + _dlclvs[u_right][f];
    }
    // D event
    values[2] = _costD + _dlclvs[u_left][e] + _dlclvs[u_right][e];
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _costS + _costL + _dlclvs[gid][f];
    values[4] = _costS + _costL + _dlclvs[gid][g];
  }
  proba = *std::max_element(values.begin(), values.end());
  if (event) {
    auto maxValueIndex = static_cast<int>(std::distance(values.begin(),
        std::max_element(values.begin(), values.end())
        ));
    if (maxValueIndex == -1 || values[maxValueIndex] == -std::numeric_limits<double>::infinity()) {
      event->type = ReconciliationEventType::EVENT_Invalid;
      return;
    }
    switch(maxValueIndex) {
    case 0:
      event->type = ReconciliationEventType::EVENT_S;
      event->cross = false;
      break;
    case 1:
      event->type = ReconciliationEventType::EVENT_S;
      event->cross = true;
      break;
    case 2:
      event->type = ReconciliationEventType::EVENT_D;
      break;
    case 3:
      event->type = ReconciliationEventType::EVENT_SL;
      event->destSpeciesNode = f;
      event->pllDestSpeciesNode = this->getSpeciesLeft(speciesNode);
      break;
    case 4:
      event->type = ReconciliationEventType::EVENT_SL;
      event->destSpeciesNode = g;
      event->pllDestSpeciesNode = this->getSpeciesRight(speciesNode);
      break;
    default:
      assert(false);
    }
  }
}
  

double ParsimonyDLModel::getGeneRootLikelihood(pll_unode_t *root) const
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


void ParsimonyDLModel::computeGeneRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    auto e = speciesNode->node_index;
    computeProbability(virtualRoot, speciesNode, _dlclvs[u][e], true);
  }
}




