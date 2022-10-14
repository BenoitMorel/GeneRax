#pragma once

#include <likelihoods/reconciliation_models/GTBaseReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/Logger.hpp>
#include <algorithm>




/*
* Implement the undated model described here:
* https://github.com/ssolo/ALE/blob/master/misc/undated.pdf
* In addition, we forbid transfers to parent species
*/
template <class REAL>
class UndatedDTLModel: public GTBaseReconciliationModel<REAL> {
public:
  UndatedDTLModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMappingp, 
      const RecModelInfo &recModelInfo): 
    GTBaseReconciliationModel<REAL>(speciesTree, 
        geneSpeciesMappingp, 
        recModelInfo),
  _transferConstraint(recModelInfo.transferConstraint) {}
  UndatedDTLModel(const UndatedDTLModel &) = delete;
  UndatedDTLModel & operator = (const UndatedDTLModel &) = delete;
  UndatedDTLModel(UndatedDTLModel &&) = delete;
  UndatedDTLModel & operator = (UndatedDTLModel &&) = delete;
  virtual ~UndatedDTLModel();
  
  // overloaded from parent
  virtual void setRates(const RatesVector &rates);
  
protected:
  // overloaded from parent
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  // overloaded from parent
  virtual void updateCLV(corax_unode_t *geneNode);
  // overload from parent
  virtual void recomputeSpeciesProbabilities();
  // overloaded from parent
  virtual REAL getGeneRootLikelihood(corax_unode_t *root) const;
  // overload from parent
  virtual void computeGeneRootLikelihood(corax_unode_t *virtualRoot);
  virtual REAL getGeneRootLikelihood(corax_unode_t *root, corax_rnode_t *speciesRoot) {
    return _dtlclvs[root->node_index + this->_maxGeneId + 1]._uq[speciesRoot->node_index];
  }
  virtual REAL getLikelihoodFactor() const;
  virtual void computeProbability(corax_unode_t *geneNode, corax_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr,
      bool stochastic = false);
private:
  // model
  std::vector<double> _PD; // Duplication probability, per branch
  std::vector<double> _PL; // Loss probability, per branch
  std::vector<double> _PT; // Transfer probability, per branch
  std::vector<double> _PS; // Speciation probability, per branch
  // SPECIES
  std::vector<double> _uE; // Probability for a gene to become extinct on each brance
  TransferConstaint _transferConstraint;
  
  /**
   *  All intermediate results needed to compute the reconciliation likelihood
   *  each gene node has one DTLCLV object
   *  Each DTLCLV gene  object is a function of the DTLCLVs of the direct children genes
   */
  struct DTLCLV {
    DTLCLV():
      _survivingTransferSums(REAL())
    {}

    DTLCLV(unsigned int speciesNumber):
      _uq(speciesNumber, REAL()),
      _correctionSum(speciesNumber, REAL()),
      _survivingTransferSums(REAL())
    {
    }
    // probability of a gene node rooted at a species node
    std::vector<REAL> _uq;
    std::vector<REAL> _correctionSum;

    // sum of transfer probabilities. Can be computed only once
    // for all species, to reduce computation complexity
    REAL _survivingTransferSums;
  };

  // Current DTLCLV values
  std::vector<DTLCLV> _dtlclvs;
  std::vector<corax_rnode_s *> _orderedSpeciations;
private:
  void getBestTransfer(corax_unode_t *parentGeneNode, 
    corax_rnode_t *originSpeciesNode,
    bool isVirtualRoot,
    corax_unode_t *&transferedGene,
    corax_unode_t *&stayingGene,
    corax_rnode_t *&recievingSpecies,
    REAL &proba, 
    bool stochastic = false);
  void getBestTransferLoss(Scenario &scenario,
      corax_unode_t *parentGeneNode, 
    corax_rnode_t *originSpeciesNode,
    corax_rnode_t *&recievingSpecies,
    REAL &proba,
    bool stochastic = false);
  unsigned int getIterationsNumber() const {return 4;}    

  REAL getCorrectedTransferSum(unsigned int geneId, unsigned int speciesId) const
  {
    switch(_transferConstraint) {
    case TransferConstaint::NONE:
      return (_dtlclvs[geneId]._survivingTransferSums
          - _dtlclvs[geneId]._uq[speciesId] * (1.0 / double(this->_allSpeciesNodes.size())))
        * _PT[speciesId];
    case TransferConstaint::PARENTS:
      return (_dtlclvs[geneId]._survivingTransferSums
          - _dtlclvs[geneId]._correctionSum[speciesId])
        * _PT[speciesId];
    case TransferConstaint::SOFTDATED:
      return _dtlclvs[geneId]._correctionSum[speciesId]
        * _PT[speciesId];
    default:
      assert(false);
    }
  }

  std::vector<corax_rnode_s *> &getSpeciesNodesToUpdate() {
    return this->_speciesNodesToUpdate;
  }

  std::vector<corax_rnode_s *> &getSpeciesNodesToUpdateSafe() {
    return this->_allSpeciesNodes;
  }
};


template <class REAL>
void UndatedDTLModel<REAL>::setInitialGeneTree(PLLUnrootedTree &tree)
{
  GTBaseReconciliationModel<REAL>::setInitialGeneTree(tree);
  DTLCLV nullCLV(this->_allSpeciesNodesCount);
  _dtlclvs = std::vector<DTLCLV>(2 * (this->_maxGeneId + 1), nullCLV);
}



template <class REAL>
void UndatedDTLModel<REAL>::setRates(const RatesVector &rates)
{
  this->_geneRoot = 0;
  assert(rates.size() == 3);
  auto &dupRates = rates[0];
  auto &lossRates = rates[1];
  auto &transferRates = rates[2];
  assert(this->_allSpeciesNodesCount == dupRates.size());
  assert(this->_allSpeciesNodesCount == lossRates.size());
  assert(this->_allSpeciesNodesCount == transferRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PT = transferRates;
  _PS.resize(this->_allSpeciesNodesCount);
  for (unsigned int e = 0; e < this->_allSpeciesNodesCount; ++e) {
    if (this->_info.noDup) {
      _PD[e] = 0.0;
    }
    auto sum = _PD[e] + _PL[e] + _PT[e] + 1.0;
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] = 1.0 / sum;
  } 
  recomputeSpeciesProbabilities();
  this->invalidateAllCLVs();
  this->invalidateAllSpeciesCLVs();
}

template <class REAL>
UndatedDTLModel<REAL>::~UndatedDTLModel() { }



template <class REAL>
void UndatedDTLModel<REAL>::recomputeSpeciesProbabilities()
{
  if (_transferConstraint == TransferConstaint::SOFTDATED) {
    _orderedSpeciations = this->_speciesTree.getOrderedSpeciations();  
  }
  _uE.resize(this->_allSpeciesNodesCount);
  for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
    _uE[speciesNode->node_index] = 0.0;
  }
  std::vector<double> transferExtinctionSums(this->_allSpeciesNodes.size(), REAL());
  for (unsigned int it = 0; it < getIterationsNumber(); ++it) {
    for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
      auto e = speciesNode->node_index;
      if (it + 1 == getIterationsNumber() && !speciesNode->left) {
        _uE[e] = _uE[e] * (1.0 - this->_fm[e]) + this->_fm[e];
        continue;
      }
      double proba = _PL[e] + (_PD[e] * _uE[e] * _uE[e]) + _PT[e] * transferExtinctionSums[e] * _uE[e];
      if (this->getSpeciesLeft(speciesNode)) {
        proba +=  _uE[this->getSpeciesLeft(speciesNode)->node_index]  * _uE[this->getSpeciesRight(speciesNode)->node_index] * _PS[e];
      }
      _uE[speciesNode->node_index] = proba;
    }
    std::fill(transferExtinctionSums.begin(), transferExtinctionSums.end(), 0.0);
    auto transferExtinctionSum = 0.0;
    double N = this->_allSpeciesNodes.size();
    if (this->_transferConstraint == TransferConstaint::NONE) {
      for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
        auto e = speciesNode->node_index;
        transferExtinctionSum += _uE[e];
      }
      transferExtinctionSum /= N;
      for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
        auto e = speciesNode->node_index;
        transferExtinctionSums[e] = transferExtinctionSum;
      }
    } else if (this->_transferConstraint == TransferConstaint::PARENTS) {
      assert(false);
    } else if (this->_transferConstraint == TransferConstaint::SOFTDATED) {
      std::vector<double> softDatedSums(N, 0.0);
      double softDatedSum = 0.0;
      for (auto leaf: this->_speciesTree.getLeaves()) {
        auto e = leaf->node_index;
        softDatedSum += _uE[e];
      }
      for (auto it = this->_orderedSpeciations.rbegin(); 
          it != this->_orderedSpeciations.rend(); ++it) {
        auto e = (*it)->node_index;
        softDatedSums[e] = softDatedSum;
        softDatedSum += _uE[e];
      }
      for (auto node: this->_allSpeciesNodes) {
        auto e = node->node_index;
        auto p = node->parent ? node->parent->node_index : e;
        transferExtinctionSums[e] = softDatedSums[p];
        if (e != p) {
          transferExtinctionSums[e] = transferExtinctionSums[e] - _uE[e]; 
        }
        transferExtinctionSums[e] /= N; 
      } 
    } else {
      assert(false);
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::updateCLV(corax_unode_t *geneNode)
{
  auto gid = geneNode->node_index; 
  auto lca = this->_geneToSpeciesLCA[gid];
  auto &clv = _dtlclvs[gid];
  auto &uq = clv._uq;
  auto &correctionSum = clv._correctionSum;
  
  auto &parentsCache = this->_speciesTree.getParentsCache(lca);
  
  /*
  auto lcaRight = lca;
  auto lcaLeft = lca;
  if (geneNode->next) {
    auto geneLeft = this->getLeft(geneNode, false);
    auto geneRight = this->getRight(geneNode, false);
    lcaRight = this->_geneToSpeciesLCA[geneLeft->node_index];
    lcaLeft = this->_geneToSpeciesLCA[geneRight->node_index];
  }
  auto &ancestorsLeft = this->_speciesTree.getAncestorssCache(lcaLeft);
  auto &ancestorsRight = this->_speciesTree.getAncestorssCache(lcaRight);
  */
 
  auto N = static_cast<double>(this->_allSpeciesNodes.size());
  std::fill(uq.begin(), uq.end(), REAL());
  std::fill(correctionSum.begin(), correctionSum.end(), REAL());
  REAL sum = REAL();
  for (auto speciesNode: getSpeciesNodesToUpdateSafe()) { 
    auto e = speciesNode->node_index;
    if (parentsCache[e]) {
    //if (ancestorsLeft[e] || ancestorsRight[e]) {
      computeProbability(geneNode, 
        speciesNode, 
        uq[e]);
    }
    sum += uq[e];
    //}
  }
  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto postOrder = this->_speciesTree.getPostOrderNodes();
    for (auto it = postOrder.rbegin();
        it != postOrder.rend(); ++it)  {
      auto speciesNode = *it;
      auto e = speciesNode->node_index;
      auto parent = speciesNode;
      while (parent) {
        auto p = parent->node_index;
        correctionSum[e] += uq[p];
        parent = parent->parent;
      }
      correctionSum[e] /= N;
    }
  }
  if (_transferConstraint == TransferConstaint::SOFTDATED) {
    std::vector<REAL> softDatedSums(N, REAL());
    std::vector<double> possibleTransfers(N, REAL());
    double currentPossibleTransfers = static_cast<double>(this->_speciesTree.getLeavesNumber()) - 1.0;
    REAL softDatedSum = REAL();
    for (auto leaf: this->_speciesTree.getLeaves()) {
      auto e = leaf->node_index;
      softDatedSum += uq[e];
      possibleTransfers[e] = currentPossibleTransfers;
    }
    for (auto it = this->_orderedSpeciations.rbegin(); 
        it != this->_orderedSpeciations.rend(); ++it) {
      auto e = (*it)->node_index;
      softDatedSums[e] = softDatedSum;
      possibleTransfers[e] = currentPossibleTransfers;
      softDatedSum += uq[e];
      currentPossibleTransfers += 1.0;
    }
    for (auto node: this->_allSpeciesNodes) {
      auto e = node->node_index;
      auto p = node->parent ? node->parent->node_index : e;
      correctionSum[e] = softDatedSums[p];
      if (e != p) {
        correctionSum[e] = correctionSum[e] - uq[e]; 
      } 
      assert (possibleTransfers[p] != 0.0);
#ifdef POSSIBLETRANSFER_NORM
      correctionSum[e] /= possibleTransfers[p];

#else
      correctionSum[e] /= N; 

#endif
      
    }
  }
  sum /= N;
  clv._survivingTransferSums = sum;
}


template <class REAL>
void UndatedDTLModel<REAL>::computeGeneRootLikelihood(corax_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  _dtlclvs[u]._survivingTransferSums = REAL();
  std::fill(_dtlclvs[u]._uq.begin(), _dtlclvs[u]._uq.end(), REAL());
  /*
  auto geneLeft = this->getLeft(virtualRoot, true);
  auto geneRight = this->getRight(virtualRoot, true);
  auto lcaRight = this->_geneToSpeciesLCA[geneLeft->node_index];
  auto lcaLeft = this->_geneToSpeciesLCA[geneRight->node_index];
  auto &ancestorsLeft = this->_speciesTree.getAncestorssCache(lcaLeft);
  auto &ancestorsRight = this->_speciesTree.getAncestorssCache(lcaRight);
  */
  for (auto speciesNode: getSpeciesNodesToUpdateSafe()) {
    unsigned int e = speciesNode->node_index;
    //if (ancestorsLeft[e] || ancestorsRight[e]) {
    computeProbability(virtualRoot, speciesNode, _dtlclvs[u]._uq[e], true);
    //}
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::computeProbability(corax_unode_t *geneNode, corax_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot,
      Scenario *scenario,
      Scenario::Event *event,
      bool stochastic)
{
  
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !this->getSpeciesLeft(speciesNode);
  
  if (event) {
    event->geneNode = gid; 
    event->speciesNode = e;
    event->type = ReconciliationEventType::EVENT_None; 
  }

  if (isSpeciesLeaf and isGeneLeaf and e == this->_geneToSpecies[gid]) {
    proba = REAL(_PS[e]);
    if (event) {
      event->label = std::string(geneNode->label);
    }
    return;
  }
  typedef std::array<REAL, 7>  ValuesArray;
  ValuesArray values;
  if (event) {
    values[0] = values[1] = values[2] = values[3] = REAL();
    values[4] = values[5] = values[6] = REAL();
  }
  proba = REAL();
  
  corax_unode_t *leftGeneNode = 0;     
  corax_unode_t *rightGeneNode = 0;     
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }

  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = this->getSpeciesLeft(speciesNode)->node_index;
    g = this->getSpeciesRight(speciesNode)->node_index;
  }

  if (not isGeneLeaf) {
    // S event
    auto u_left = leftGeneNode->node_index;
    auto u_right = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      //  speciation event
      values[0] = _dtlclvs[u_left]._uq[f];
      values[1] = _dtlclvs[u_left]._uq[g];
      values[0] *= _dtlclvs[u_right]._uq[g];
      values[1] *= _dtlclvs[u_right]._uq[f];
      values[0] *= _PS[e]; 
      values[1] *= _PS[e]; 
      scale(values[0]);
      scale(values[1]);
      proba += values[0];
      proba += values[1];
    }
    // D event
    values[2] = _dtlclvs[u_left]._uq[e];
    values[2] *= _dtlclvs[u_right]._uq[e];
    values[2] *= _PD[e];
    scale(values[2]);
    proba += values[2];
    
    // T event
    values[5] = getCorrectedTransferSum(u_left, e);
    values[5] *= _dtlclvs[u_right]._uq[e];
    scale(values[5]);
    values[6] = getCorrectedTransferSum(u_right, e);
    values[6] *= _dtlclvs[u_left]._uq[e];
    scale(values[6]);
    proba += values[5];
    proba += values[6];
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _dtlclvs[gid]._uq[f];
    values[3] *= (_uE[g] * _PS[e]);
    scale(values[3]);
    values[4] = _dtlclvs[gid]._uq[g];
    values[4]*= _uE[f] * _PS[e];
    scale(values[4]);
    proba += values[3];
    proba += values[4];
  }
  if (event) {
    assert(scenario);
    corax_unode_t *transferedGene = 0;
    corax_unode_t *stayingGene = 0;
    corax_rnode_t *recievingSpecies = 0;
    values[5] = values[6] = REAL(); // invalidate these ones
    if (!isGeneLeaf) {
      getBestTransfer(geneNode, speciesNode, isVirtualRoot, 
          transferedGene, stayingGene, recievingSpecies, values[5], stochastic);
    }
    int maxValueIndex = 0;
    if (!stochastic) {
      maxValueIndex =static_cast<unsigned int>(std::distance(values.begin(),
          std::max_element(values.begin(), values.end())
          ));
      
      assert(values[maxValueIndex] != REAL());
    } else {
      maxValueIndex = sampleIndex<ValuesArray, REAL>(values);
    }
    if (-1 == maxValueIndex || values[maxValueIndex] == REAL()) {
      event->type = ReconciliationEventType::EVENT_Invalid;
      return;
    }
    switch(maxValueIndex) {
    case 0:
      event->type = ReconciliationEventType::EVENT_S;
      event->leftGeneIndex = leftGeneNode->node_index;
      event->rightGeneIndex = rightGeneNode->node_index;
      break;
    case 1:
      event->type = ReconciliationEventType::EVENT_S;
      event->leftGeneIndex = rightGeneNode->node_index;
      event->rightGeneIndex = leftGeneNode->node_index;
      break;
    case 2:
      event->type = ReconciliationEventType::EVENT_D;
      event->leftGeneIndex = leftGeneNode->node_index;
      event->rightGeneIndex = rightGeneNode->node_index;
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
      //assert (values[5] != proba);
      event->type = ReconciliationEventType::EVENT_T;
      event->leftGeneIndex = stayingGene->node_index;
      event->rightGeneIndex = transferedGene->node_index;
      event->destSpeciesNode = recievingSpecies->node_index;
      event->pllTransferedGeneNode = transferedGene;
      event->pllDestSpeciesNode = recievingSpecies;
      break;
    case 6:
      assert(false);
      break;
    default:
      assert(false);
    };
  }
}



template <class REAL>
REAL UndatedDTLModel<REAL>::getGeneRootLikelihood(corax_unode_t *root) const
{
  REAL sum = REAL();
  auto u = root->node_index + this->_maxGeneId + 1;
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    sum += _dtlclvs[u]._uq[e];
  }
  return sum;
}

template <class REAL>
REAL UndatedDTLModel<REAL>::getLikelihoodFactor() const
{
  REAL factor(0.0);
  for (auto speciesNode: this->_allSpeciesNodes) {
    auto e = speciesNode->node_index;
    factor += (REAL(1.0) - REAL(_uE[e]));
  }
  return factor;
      
}

template <class REAL>
void UndatedDTLModel<REAL>::getBestTransfer(corax_unode_t *parentGeneNode, 
  corax_rnode_t *originSpeciesNode,
  bool isVirtualRoot,
  corax_unode_t *&transferedGene,
  corax_unode_t *&stayingGene,
  corax_rnode_t *&recievingSpecies,
  REAL &proba, 
  bool stochastic)
{
  unsigned int speciesNumber = this->_speciesTree.getNodesNumber();;
  proba = REAL();
  auto e = originSpeciesNode->node_index;
  std::unordered_set<unsigned int> parents;
  if (_transferConstraint == TransferConstaint::PARENTS) {
    auto parent = originSpeciesNode;
    while (parent) {
      parents.insert(parent->node_index);
      parent = parent->parent;
    }
  }
  auto u_left = this->getLeft(parentGeneNode, isVirtualRoot);
  auto u_right = this->getRight(parentGeneNode, isVirtualRoot);
  std::vector<REAL> transferProbas(speciesNumber * 2, REAL());
  double factor = _PT[e] / static_cast<double>(speciesNumber);
  for (auto species: this->_allSpeciesNodes) {
    auto h = species->node_index;
    if (_transferConstraint == TransferConstaint::PARENTS) {
      if (parents.end() != parents.find(h)) {
        continue;
      }
    }
    if (_transferConstraint == TransferConstaint::NONE) {
      if (h == e) {
        continue;
      }
    }

    transferProbas[h] = (_dtlclvs[u_left->node_index]._uq[h] 
        * _dtlclvs[u_right->node_index]._uq[e]) * factor;
    transferProbas[h + speciesNumber] = (_dtlclvs[u_right->node_index]._uq[h] 
        * _dtlclvs[u_left->node_index]._uq[e]) * factor;
  }
  if (stochastic) {
    // stochastic sample: proba will be set to the sum of probabilities
    proba = REAL();
    for (auto &value: transferProbas) {
      proba += value;
    }
    auto bestIndex = sampleIndex<std::vector<REAL>, REAL>(transferProbas);
    if (bestIndex == -1) {
      proba = REAL();
      return;
    }
    bool left = static_cast<unsigned int>(bestIndex) < speciesNumber;
    transferedGene = left ? u_left : u_right;
    stayingGene = !left ? u_left : u_right;
    // I am not sure I can find the species from
    // its index in the species array, so I keep it safe here
    for (auto species: this->_allSpeciesNodes) {
      if (species->node_index == bestIndex % speciesNumber) {
        recievingSpecies = species;
      }
    }
  } else {
    // find the max
    REAL sum = REAL();
    for (auto species: this->_allSpeciesNodes) {
      auto h = species->node_index;
      if (parents.end() != parents.find(h)) {
        continue;
      }
      sum += transferProbas[h] + transferProbas[h + speciesNumber];
      if (proba < transferProbas[h]) {
        proba = transferProbas[h];
        transferedGene = u_left;
        stayingGene = u_right;
        recievingSpecies = species;
      }
      if (proba < transferProbas[h + speciesNumber]) {
        proba = transferProbas[h + speciesNumber];
        transferedGene = u_right;
        stayingGene = u_left;
        recievingSpecies = species;
      }
    }
  }
}

template <class REAL>
void UndatedDTLModel<REAL>::getBestTransferLoss(Scenario &scenario,
   corax_unode_t *parentGeneNode, 
  corax_rnode_t *originSpeciesNode,
  corax_rnode_t *&recievingSpecies,
  REAL &proba, 
  bool stochastic)
{
  proba = REAL();
  auto e = originSpeciesNode->node_index;
  auto u = parentGeneNode->node_index;
  
  unsigned int speciesNumber = this->_speciesTree.getNodesNumber();;
  std::vector<REAL> transferProbas(speciesNumber, REAL());
  REAL factor = _uE[e] * (_PT[e] / static_cast<double>(this->_allSpeciesNodes.size()));
  for (auto species: this->_allSpeciesNodes) {
    auto h = species->node_index;
    if (h == e) {
      continue;
    }
    transferProbas[h] = _dtlclvs[u]._uq[h] * factor;
  }
  if (!stochastic) {
    for (auto species: this->_allSpeciesNodes) {
      auto h = species->node_index;
      if (proba < transferProbas[h]) {
        if (!scenario.isBlacklisted(u, h)) {
          scenario.blackList(u, h);
          proba = transferProbas[h];
          recievingSpecies = species;
        }
      }
    } 
  } else {
    proba = REAL();
    for (auto p: transferProbas) {
      proba += p;
    }
    unsigned int h = 0;
    do {

      int bestIndex = sampleIndex<std::vector<REAL>, REAL>(transferProbas);
      if (bestIndex == -1) {
        proba = REAL();
        return;
      }
      for (auto species: this->_allSpeciesNodes) {
        h = species->node_index;
        if (static_cast<unsigned int>(bestIndex) == h) {
          recievingSpecies = species;
          transferProbas[h] = REAL(); // in case it's blacklisted, avoid infinite loop
          break;
        }
      }
    } while (scenario.isBlacklisted(u, h));
    scenario.blackList(u, h);
  }
}
