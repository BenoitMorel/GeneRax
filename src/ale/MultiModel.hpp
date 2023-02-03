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

struct Highway {
  Highway():src(nullptr), dest(nullptr), proba(0.0) {}

  Highway(corax_rnode_t *src, corax_rnode_t *dest):
    src(src),
    dest(dest),
    proba(0.1)
  {}
  corax_rnode_t *src;
  corax_rnode_t *dest;
  double proba;
  
  friend std::ostream& operator<<(std::ostream& os, const Highway &h) {
    os << "(" << h.src->label << " -> " << h.dest->label << ")";
    return os;
  }
};
  


class MultiModel: public BaseReconciliationModel {
public:
  MultiModel(PLLRootedTree &speciesTree,
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &ccpFile):
    BaseReconciliationModel(speciesTree,
        geneSpeciesMapping,
        info)
  {
    _ccp.unserialize(ccpFile);
    mapGenesToSpecies();
  }
  virtual ~MultiModel() {}

  const ConditionalClades &getCCP() const {return _ccp;}

  virtual void setHighways(const std::vector<Highway> &highways) {(void)(highways);}

  virtual void onSpeciesDatesChange() {}
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
        this->getAllSpeciesNodeNumber(), 0);
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
      const std::string &ccpFile):
    MultiModel(speciesTree,
        geneSpeciesMapping,
        info,
        ccpFile){}
  virtual ~MultiModelTemplate() {}
  virtual bool inferMLScenario(Scenario &scenario, 
    bool stochastic = false);
protected:
  virtual corax_rnode_t *sampleSpeciesNode(unsigned int &category) = 0;

private:
  virtual void computeProbability(CID cid, 
    corax_rnode_t *speciesNode, 
    size_t category,
    REAL &proba,
    ReconciliationCell<REAL> *recCell = nullptr) = 0;
  bool backtrace(unsigned int cid, 
      corax_rnode_t *speciesRoot,
      corax_unode_t *geneNode,
      unsigned int category,
      Scenario &scenario,
      bool stochastic); 
};




template <class REAL>
bool MultiModelTemplate<REAL>::inferMLScenario(Scenario &scenario, 
    bool stochastic) {
  assert(stochastic); 
  unsigned int category = 0;
  auto rootCID = this->_ccp.getCladesNumber() - 1;
  auto speciesRoot = sampleSpeciesNode(category);
  auto rootIndex = 2 * _ccp.getLeafNumber(); 
  scenario.setVirtualRootIndex(rootIndex);
  scenario.setSpeciesTree(&this->_speciesTree);
  auto virtualGeneRoot = scenario.generateVirtualGeneRoot();
  scenario.setGeneRoot(virtualGeneRoot);
  auto res = backtrace(rootCID, 
      speciesRoot,
      virtualGeneRoot,
      category,
      scenario, 
      stochastic);

  return res;
}


template <class REAL>
bool MultiModelTemplate<REAL>::backtrace(unsigned int cid, 
    corax_rnode_t *speciesNode,
    corax_unode_t *geneNode,
    unsigned int category,
    Scenario &scenario,
    bool stochastic)
{
  REAL proba;
  computeProbability(cid, speciesNode, category, proba);
  ReconciliationCell<REAL> recCell;
  recCell.maxProba = proba * Random::getProba();
  computeProbability(cid, speciesNode, category, proba, &recCell);
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
  auto c = category;

  switch(recCell.event.type) {
  case ReconciliationEventType::EVENT_S:
    ok &= backtrace(leftCid, this->getSpeciesLeft(speciesNode), leftGeneNode, c, scenario, stochastic); 
    ok &= backtrace(rightCid, this->getSpeciesRight(speciesNode), rightGeneNode, c, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_D:
    ok &= backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic); 
    ok &= backtrace(rightCid, speciesNode, rightGeneNode, c, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_SL:
    ok &= backtrace(cid, recCell.event.pllDestSpeciesNode, geneNode, c, scenario, stochastic); 
    break;
  case ReconciliationEventType::EVENT_T:
    // source species
    ok &= backtrace(leftCid, speciesNode, leftGeneNode, c, scenario, stochastic);
    // dest species
    ok &= backtrace(rightCid, recCell.event.pllDestSpeciesNode, rightGeneNode, c, scenario, stochastic);
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






