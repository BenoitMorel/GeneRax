#pragma once

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <util/Scenario.hpp>
#include <cmath>



static const double WORST_COST = -std::numeric_limits<double>::infinity();

class ParsimonyDTLModel: public AbstractReconciliationModel<double> {
public:
  ParsimonyDTLModel(PLLRootedTree &speciesTree, const GeneSpeciesMapping &geneSpeciesMappingp, 
      bool rootedGeneTree,
      double minGeneBranchLength,
      bool pruneSpeciesTree):
    AbstractReconciliationModel<double>(speciesTree, geneSpeciesMappingp, rootedGeneTree, minGeneBranchLength, pruneSpeciesTree),
    _costS(-0.0001),
    _costD(-1.0),
    _costL(-0.0001),
    _costT(-5.0)
    {}
  
  
  ParsimonyDTLModel(const ParsimonyDTLModel &) = delete;
  ParsimonyDTLModel & operator = (const ParsimonyDTLModel &) = delete;
  ParsimonyDTLModel(ParsimonyDTLModel &&) = delete;
  ParsimonyDTLModel & operator = (ParsimonyDTLModel &&) = delete;
  virtual ~ParsimonyDTLModel() {}
  
  // overloaded from parent
  virtual void setRates(const RatesVector &){}
  
  unsigned int getIterationsNumber() const {return 2;}    
protected:
  // overload from parent
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  // overload from parent
  virtual void updateCLV(pll_unode_t *geneNode);
  // overload from parent
  virtual double getGeneRootLikelihood(pll_unode_t *root) const;
  virtual double  getGeneRootLikelihood(pll_unode_t *root, 
      pll_rnode_t *speciesRoot) {
    return _dtlclvs[root->node_index + this->_maxGeneId + 1]
      .cost[speciesRoot->node_index];
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
  double _costT;
  
  // uq[geneId][speciesId] = parsimony cost
  // of a gene node rooted at a species node
  // to produce the subtree of this gene node
  
  struct DTLCLV {
    std::vector<double> cost;
    double bestCost;
    unsigned int bestSpecies;
    DTLCLV():
      cost(),
      bestCost(WORST_COST),
      bestSpecies(-1){}

    DTLCLV(unsigned int speciesNumber):
      cost(speciesNumber, WORST_COST),
      bestCost(WORST_COST),
      bestSpecies(-1) {}

    void reset() {
      std::fill(cost.begin(), cost.end(), WORST_COST);
      bestCost = WORST_COST;
      bestSpecies = -1;
    }

  };
  std::vector<DTLCLV> _dtlclvs;
 
private:
  std::vector<pll_rnode_s *> &getSpeciesNodesToUpdate() {
    return this->_allSpeciesNodes;
  }

};



void ParsimonyDTLModel::setInitialGeneTree(PLLUnrootedTree &tree)
{
  AbstractReconciliationModel<double>::setInitialGeneTree(tree);
  assert(this->_allSpeciesNodesCount);
  assert(this->_maxGeneId);
  DTLCLV initialCLV(this->_allSpeciesNodesCount);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_maxGeneId + 1), initialCLV);

}


void ParsimonyDTLModel::updateCLV(pll_unode_t *geneNode)
{
  assert(geneNode);
  auto gid = geneNode->node_index;
  _dtlclvs[gid].reset();
  for (unsigned int i = 0; i < getIterationsNumber(); ++i) {
    for (auto speciesNode: getSpeciesNodesToUpdate()) {
      computeProbability(geneNode, 
          speciesNode, 
          _dtlclvs[gid].cost[speciesNode->node_index]);
    }
    auto max = std::max_element(_dtlclvs[gid].cost.begin(),
        _dtlclvs[gid].cost.end());

    _dtlclvs[gid].bestCost = *max;
    _dtlclvs[gid].bestSpecies = std::distance(_dtlclvs[gid].cost.begin(), max);
  }
}



void ParsimonyDTLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
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


  proba = WORST_COST;
  if (isSpeciesLeaf and isGeneLeaf) {
    // present
    if (e == this->_geneToSpecies[gid]) {
      proba = 0.0;
    }
    return;
  }
  typedef std::array<double, 8>  ValuesArray;
  ValuesArray values;
  for (unsigned int i = 0; i < 8; ++i) {
    values[i] = WORST_COST;
  }
  unsigned int u_left = -1; 
  unsigned int u_right = -1; 
  if (not isGeneLeaf) {
    u_left = leftGeneNode->node_index;
    u_right = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      values[0] = _costS
        + _dtlclvs[u_left].cost[f] 
        + _dtlclvs[u_right].cost[g];
      values[1] = _costS 
        + _dtlclvs[u_left].cost[g] 
        + _dtlclvs[u_right].cost[f];
    }
    // D event
    values[2] = _costD
      + _dtlclvs[u_left].cost[e] 
      + _dtlclvs[u_right].cost[e];
    // T event
    values[5] = _costT
      + _dtlclvs[u_left].bestCost
      + _dtlclvs[u_right].cost[e];
    values[6] = _costT
      + _dtlclvs[u_right].bestCost
      + _dtlclvs[u_left].cost[e];
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _costS + _costL + _dtlclvs[gid].cost[f];
    values[4] = _costS + _costL + _dtlclvs[gid].cost[g];
  }
  // TL event
  values[7] = _costT + _costL + _dtlclvs[gid].bestCost;
  proba = *std::max_element(values.begin(), values.end());
  /*
  auto maxValueIndex = static_cast<int>(std::distance(values.begin(),
        std::max_element(values.begin(), values.end())));
  Logger::info << maxValueIndex << std::endl; 
        */
  if (event) {
    auto maxValueIndex = static_cast<int>(std::distance(values.begin(),
        std::max_element(values.begin(), values.end())
        ));
    if (maxValueIndex == -1 || values[maxValueIndex] == WORST_COST) {
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
    case 5:
      event->type = ReconciliationEventType::EVENT_T;
      event->transferedGeneNode = u_left;
      event->destSpeciesNode = _dtlclvs[u_left].bestSpecies;
      event->pllTransferedGeneNode = leftGeneNode;
      event->pllDestSpeciesNode = this->_allSpeciesNodes[event->destSpeciesNode];
      break;
    case 6:
      event->type = ReconciliationEventType::EVENT_T;
      event->transferedGeneNode = u_right;
      event->destSpeciesNode = _dtlclvs[u_right].bestSpecies;
      event->pllTransferedGeneNode = rightGeneNode;
      event->pllDestSpeciesNode = this->_allSpeciesNodes[event->destSpeciesNode];
      break;
    case 7:
      event->type = ReconciliationEventType::EVENT_TL;
      event->transferedGeneNode = gid;
      event->destSpeciesNode = _dtlclvs[gid].bestSpecies;
      event->pllTransferedGeneNode = geneNode;
      event->pllDestSpeciesNode = this->_allSpeciesNodes[event->destSpeciesNode];
      break;
    default:
      assert(false);
    }
  }
}
  

double ParsimonyDTLModel::getGeneRootLikelihood(pll_unode_t *root) const
{
  double max =  WORST_COST;
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    if (max < _dtlclvs[u].cost[e]) {
      max = _dtlclvs[u].cost[e];
    }
  }
  return max;
}


void ParsimonyDTLModel::computeGeneRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  for (auto speciesNode: getSpeciesNodesToUpdate()) {
    auto e = speciesNode->node_index;
    computeProbability(virtualRoot, speciesNode, _dtlclvs[u].cost[e], true);
  }
}




