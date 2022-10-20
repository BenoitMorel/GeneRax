#pragma once

#include <ccp/ConditionalClades.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>
#include <IO/LibpllParsers.hpp>


class GeneSpeciesMapping;
class PLLRootedTree;
 
template<class REAL>
struct ReconciliationCell {
  Scenario::Event event;
  REAL maxProba;
};


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

  const ConditionalClades &getCCP() const {return _ccp;}
protected:
  void mapGenesToSpecies() {
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
  ConditionalClades _ccp;
};

template<class REAL>
class MultiModelTemplate: public MultiModel {
public:
  MultiModelTemplate(PLLRootedTree &speciesTree,
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile):
    MultiModel(speciesTree,
        geneSpeciesMapping,
        info,
        geneTreesFile){}
  virtual ~MultiModelTemplate() {}
  virtual bool inferMLScenario(Scenario &scenario, 
    bool stochastic = false);
protected:
  virtual corax_rnode_t *sampleSpeciesNode() = 0;

private:
  virtual void computeProbability(CID cid, 
    corax_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell<REAL> *recCell = nullptr) = 0;
  bool backtrace(unsigned int cid, 
      corax_rnode_t *speciesRoot,
      corax_unode_t *geneNode,
      Scenario &scenario,
      bool stochastic); 
};




template <class REAL>
bool MultiModelTemplate<REAL>::inferMLScenario(Scenario &scenario, 
    bool stochastic) {
   
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto speciesRoot = sampleSpeciesNode();
  auto rootIndex = 2 * _ccp.getLeafNumber(); 
  scenario.setVirtualRootIndex(rootIndex);
  scenario.setSpeciesTree(&this->_speciesTree);
  auto virtualGeneRoot = scenario.generateVirtualGeneRoot();
  scenario.setGeneRoot(virtualGeneRoot);
  auto res = backtrace(rootCID, 
      speciesRoot,
      virtualGeneRoot, 
      scenario, 
      stochastic);

  return res;
}


template <class REAL>
bool MultiModelTemplate<REAL>::backtrace(unsigned int cid, 
    corax_rnode_t *speciesNode,
    corax_unode_t *geneNode,
    Scenario &scenario,
    bool stochastic)
{
  REAL proba;
  computeProbability(cid, speciesNode, proba);
  ReconciliationCell<REAL> recCell;
  recCell.maxProba = proba * Random::getProba();
  computeProbability(cid, speciesNode, proba, &recCell);
  if (scenario.getGeneNodeBuffer().size() == 1) {
    recCell.event.geneNode = scenario.getVirtualRootIndex();
  } else {
    recCell.event.geneNode = geneNode->node_index;
  }
  
  corax_unode_t *leftGeneNode;
  corax_unode_t *rightGeneNode;
  auto leftCid = recCell.event.leftGeneIndex;
  auto rightCid = recCell.event.rightGeneIndex;
  if (recCell.event.type == ReconciliationEventType::EVENT_S 
      || recCell.event.type == ReconciliationEventType::EVENT_D
      || recCell.event.type == ReconciliationEventType::EVENT_T) {
    scenario.generateGeneChildren(geneNode, leftGeneNode, rightGeneNode);
    recCell.event.leftGeneIndex = leftGeneNode->node_index;
    recCell.event.rightGeneIndex = rightGeneNode->node_index;
  }
  scenario.addEvent(recCell.event);

  bool ok = true;
  std::string label;

  switch(recCell.event.type) {
  case ReconciliationEventType::EVENT_S:
    ok &= backtrace(leftCid, this->getSpeciesLeft(speciesNode), leftGeneNode, scenario, stochastic); 
    ok &= backtrace(rightCid, this->getSpeciesRight(speciesNode), rightGeneNode, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_D:
    ok &= backtrace(leftCid, speciesNode, leftGeneNode, scenario, stochastic); 
    ok &= backtrace(rightCid, speciesNode, rightGeneNode, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_SL:
    ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, geneNode, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_T:
    // source species
    ok &= backtrace(leftCid, speciesNode, leftGeneNode, scenario, stochastic);
    // dest species
    ok &= backtrace(rightCid, recCell.event.pllDestSpeciesNode, rightGeneNode, scenario, stochastic);
    break;
  case ReconciliationEventType::EVENT_TL:
    assert(false);
    break;
  case ReconciliationEventType::EVENT_None:
    label = _ccp.getCidToLeaves().at(cid);
    geneNode->label = new char[label.size() + 1];
    memcpy(geneNode->label, label.c_str(), label.size() + 1);
    break;
  default:
    ok  = false;
  }
  return ok;
}






