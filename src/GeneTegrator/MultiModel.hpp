#pragma once

#include <ccp/ConditionalClades.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>


class GeneSpeciesMapping;
class PLLRootedTree;
 
template<class REAL>
struct ReconciliationCell {
  Scenario::Event event;
  REAL maxProba;
};

template<class REAL>
class MultiModel: public BaseReconciliationModel {
public:
  MultiModel(PLLRootedTree &speciesTree,
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile):
    BaseReconciliationModel(speciesTree,
        geneSpeciesMapping,
        info),
    _ccp(geneTreesFile)
  {
    mapGenesToSpecies();
  }
  virtual ~MultiModel() {}
  virtual bool inferMLScenario(Scenario &scenario, 
    bool stochastic = false);

protected:
  void mapGenesToSpecies();
  virtual void computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell<REAL> *recCell = nullptr) = 0;
  bool backtrace(unsigned int cid, 
      pll_rnode_t *speciesRoot,
      Scenario &scenario,
      bool stochastic); 



  ConditionalClades _ccp;
};

template <class REAL>
void MultiModel<REAL>::mapGenesToSpecies()
{
  const auto &cidToLeaves = _ccp.getCidToLeaves();
  this->_speciesNameToId.clear();
  this->_geneToSpecies.resize(_ccp.getCladesNumber());
  for (auto node: this->_allSpeciesNodes) {
    if (!node->left) {
      this->_speciesNameToId[node->label] = node->node_index;
    }
  }
  this->_speciesCoverage = std::vector<unsigned int>(
      this->_allSpeciesNodesCount, 0);
  for (auto p: cidToLeaves) {
    auto cid = p.first;
    const auto &geneName = cidToLeaves.at(cid);
    const auto &speciesName = this->_geneNameToSpeciesName[geneName];
    this->_geneToSpecies[cid] = this->_speciesNameToId[speciesName];
    this->_speciesCoverage[this->_geneToSpecies[cid]]++;
  }
}

template <class REAL>
bool MultiModel<REAL>::inferMLScenario(Scenario &scenario, 
    bool stochastic) {
  // TODO: sample the species root!
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto speciesRoot = this->_speciesTree.getRoot();
  return backtrace(rootCID, 
      speciesRoot,
      scenario, stochastic);
}

static double dRand()
{
  return  (double)rand() / RAND_MAX;
}

template <class REAL>
bool MultiModel<REAL>::backtrace(unsigned int cid, 
    pll_rnode_t *speciesNode,
    Scenario &scenario,
    bool stochastic)
{
  REAL proba;
  computeProbability(cid, speciesNode, proba);
  ReconciliationCell<REAL> recCell;
  recCell.maxProba = proba * dRand();
  computeProbability(cid, speciesNode, proba, &recCell);
  scenario.addEvent(recCell.event);
  bool ok = true;
  switch(recCell.event.type) {
  case ReconciliationEventType::EVENT_S:
    ok &= backtrace(recCell.event.leftGeneIndex, this->getSpeciesLeft(speciesNode), scenario, stochastic); 
    ok &= backtrace(recCell.event.rightGeneIndex, this->getSpeciesRight(speciesNode), scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_D:
    ok &= backtrace(recCell.event.leftGeneIndex, speciesNode, scenario, stochastic); 
    ok &= backtrace(recCell.event.rightGeneIndex, speciesNode, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_SL:
    ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_T:
    // source species
    ok &= backtrace(recCell.event.leftGeneIndex, speciesNode, scenario, stochastic);
    // dest species
    ok &= backtrace(recCell.event.rightGeneIndex, recCell.event.pllDestSpeciesNode, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_TL:
    assert(false);
    break;
  case ReconciliationEventType::EVENT_None:
    break;
  default:
    ok  = false;
  }
  return ok;
}






