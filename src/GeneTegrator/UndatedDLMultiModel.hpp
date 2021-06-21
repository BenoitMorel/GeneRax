#pragma once

#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include <ccp/ConditionalClades.hpp>
#include <maths/ScaledValue.hpp>
#include <likelihoods/reconciliation_models/BaseReconciliationModel.hpp>

class RecModelInfo;
double log(ScaledValue v);

template <class REAL>
class UndatedDLMultiModel: public BaseReconciliationModel {
public: 
  UndatedDLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &info,
      const std::string &geneTreesFile);

  virtual ~UndatedDLMultiModel() {}

  virtual void setRates(const RatesVector &) {};
  virtual double computeLogLikelihood();
  virtual bool inferMLScenario(Scenario &, bool stochastic = false) {
    (void)stochastic;
    return false;}
  
private:
  
  
  ConditionalClades _ccp;
  std::vector<double> _PD; // Duplication probability, per species branch
  std::vector<double> _PL; // Loss probability, per species branch
  std::vector<double> _PS; // Speciation probability, per species branch
  std::vector<double> _uE; // Extinction probability, per species branch
  using DLCLV = std::vector<REAL>;
  std::vector<DLCLV> _dlclvs;

  struct ReconciliationCell {
    Scenario::Event event;
    REAL maxProba;
  };

  REAL getLikelihoodFactor() const;
  virtual void recomputeSpeciesProbabilities();
  void computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell *recCell = nullptr);

  void mapGenesToSpecies();
  
};

template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &info,
    const std::string &geneTreesFile):
  BaseReconciliationModel(speciesTree,
      geneSpeciesMapping,
      info),
  _ccp(geneTreesFile),
  _PD(speciesTree.getNodesNumber(), 0.2),
  _PL(speciesTree.getNodesNumber(), 0.2),
  _PS(speciesTree.getNodesNumber(), 1.0),
  _uE(speciesTree.getNodesNumber(), 0.0)
{
  std::vector<REAL> zeros(speciesTree.getNodesNumber(), REAL());
  _dlclvs = std::vector<std::vector<REAL> >(
      _ccp.getCladesNumber(), zeros);
  for (unsigned int e = 0; e < _speciesTree.getNodesNumber(); ++e) {
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  mapGenesToSpecies();
}

template <class REAL>
double UndatedDLMultiModel<REAL>::computeLogLikelihood()
{ 
  if (_ccp.skip()) {
    return 0.0;
  }
  beforeComputeLogLikelihood();
  std::vector<REAL> zeros(_speciesTree.getNodesNumber(), REAL());
  _dlclvs = std::vector<std::vector<REAL> >(
      _ccp.getCladesNumber(), zeros);
  for (CID cid = 0; cid < _ccp.getCladesNumber(); ++cid) {
    for (auto speciesNode: _allSpeciesNodes) {
      computeProbability(cid, 
          speciesNode, 
          _dlclvs[cid][speciesNode->node_index]);
    }
  }
  auto rootCID = _ccp.getCladesNumber() - 1;
  REAL res = REAL();
  for (auto speciesNode: _allSpeciesNodes) {
    res += _dlclvs[rootCID][speciesNode->node_index];
  }
  // the root correction makes sure that UndatedDLMultiModel and
  // UndatedDL model are equivalent when there is one tree per
  // family: the UndatedDLMultiModel integrates over all possible
  // roots and adds a 1/numberOfGeneRoots weight that is not
  // present un the UndatedDL, so we multiply back here
  REAL rootCorrection(double(_ccp.getRootsNumber()));
  return log(res) - log(getLikelihoodFactor()) + log(rootCorrection);
}

template <class REAL>
void UndatedDLMultiModel<REAL>::recomputeSpeciesProbabilities()
{
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (getSpeciesLeft(speciesNode)) {
      c += _PS[e] * _uE[getSpeciesLeft(speciesNode)->node_index]  * 
        _uE[getSpeciesRight(speciesNode)->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    _uE[speciesNode->node_index] = proba;
  }
}

template <class REAL>
void UndatedDLMultiModel<REAL>::computeProbability(CID cid, 
    pll_rnode_t *speciesNode, 
    REAL &proba,
    ReconciliationCell *recCell
    )
{
  proba = REAL();
  bool isSpeciesLeaf = !getSpeciesLeft(speciesNode);
  auto e = speciesNode->node_index;
  REAL maxProba = REAL();
  if (recCell) {
    recCell->event.geneNode = cid; 
    recCell->event.speciesNode = e;
    recCell->event.type = ReconciliationEventType::EVENT_None; 
    maxProba = recCell->maxProba;
  }
  // terminal gene and species nodes
  if (_ccp.isLeaf(cid) && isSpeciesLeaf) {
    if (_geneToSpecies[cid] == e) {
      proba = REAL(_PS[e]);
    }
    return;
  }
  REAL temp;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = getSpeciesLeft(speciesNode)->node_index;
    g = getSpeciesRight(speciesNode)->node_index;
  }
  
  for (const auto &cladeSplit: _ccp.getCladeSplits(cid)) {
    auto cidLeft = cladeSplit.left; 
    auto cidRight = cladeSplit.right;
    auto freq = cladeSplit.frequency;
    if (not isSpeciesLeaf) {
      // S event;
      temp = _dlclvs[cidLeft][f] * _dlclvs[cidRight][g] * (_PS[e] * freq); 
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidLeft;
        recCell->event.rightGeneIndex = cidRight;
        return;
      }
      temp = _dlclvs[cidRight][f] * _dlclvs[cidLeft][g] * (_PS[e] * freq); 
      scale(temp);
      proba += temp;
      if (recCell && proba > maxProba) {
        recCell->event.type = ReconciliationEventType::EVENT_S;
        recCell->event.leftGeneIndex = cidRight;
        recCell->event.rightGeneIndex = cidLeft;
        return;
      }
    }
    temp = _dlclvs[cidLeft][e] * _dlclvs[cidRight][e] * (_PD[e] * freq);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_D;
      recCell->event.leftGeneIndex = cidLeft;
      recCell->event.rightGeneIndex = cidRight;
      return;
    }
  }
  if (not isSpeciesLeaf) {
    // SL event
    temp = _dlclvs[cid][f] * (_uE[g] * _PS[e]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = f;
      return;
    }
    
    temp = _dlclvs[cid][g] * (_uE[f] * _PS[e]);
    scale(temp);
    proba += temp;
    if (recCell && proba > maxProba) {
      recCell->event.type = ReconciliationEventType::EVENT_SL;
      recCell->event.destSpeciesNode = g;
      return;
    }
  }
  // DL event
  //proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
  if (recCell) {
    // we haven't sampled any event...
    assert(false);
  }
}
  
template <class REAL>
REAL UndatedDLMultiModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
}


  template <class REAL>
void UndatedDLMultiModel<REAL>::mapGenesToSpecies()
{
  const auto &cidToLeaves = _ccp.getCidToLeaves();
  _speciesNameToId.clear();
  this->_geneToSpecies.resize(_ccp.getCladesNumber());
  for (auto node: _allSpeciesNodes) {
    if (!node->left) {
      _speciesNameToId[node->label] = node->node_index;
    }
  }
  this->_speciesCoverage = std::vector<unsigned int>(
      this->_allSpeciesNodesCount, 0);
  for (auto p: cidToLeaves) {
    auto cid = p.first;
    const auto &geneName = cidToLeaves.at(cid);
    const auto &speciesName = _geneNameToSpeciesName[geneName];
    _geneToSpecies[cid] = _speciesNameToId[speciesName];
    _speciesCoverage[_geneToSpecies[cid]]++;
  }
}

