#include "UndatedDTLModel.hpp"
#include <IO/Logger.hpp>
#include <algorithm>


const unsigned int CACHE_SIZE = 100000;
const unsigned int IT = 5;

UndatedDTLModel::UndatedDTLModel()
{
  _maxGeneId = 1;
}

#define ASSERT_PROBA(x) assert(x.isProba());

void UndatedDTLModel::setInitialGeneTree(pll_utree_t *tree)
{
  AbstractReconciliationModel::setInitialGeneTree(tree);
  std::vector<ScaledValue> zeros(speciesNodesCount_);
  _uq = std::vector<std::vector<ScaledValue> >(2 * (_maxGeneId + 1),zeros);
  _survivingTransferSums = std::vector<ScaledValue>(2 * (_maxGeneId + 1));
  _ancestralCorrection = std::vector<std::vector<ScaledValue> >(2 * (_maxGeneId + 1),zeros);
}

void UndatedDTLModel::updateTransferSums(ScaledValue &transferSum,
    std::vector<ScaledValue> &ancestralCorrection,
    const std::vector<ScaledValue> &probabilities)
{
  transferSum = ScaledValue();
  for (int i = static_cast<int>(speciesNodes_.size()) - 1; i >= 0; --i) {
    auto speciesNode = speciesNodes_[static_cast<unsigned int>(i)];
    auto e = speciesNode->node_index;
    ancestralCorrection[e] = probabilities[e] * _PT[e];
    if (speciesNode->parent) {
      auto p = speciesNode->parent->node_index;
      ancestralCorrection[e] += ancestralCorrection[p];
    }
  }
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    ancestralCorrection[e] /= double(speciesNodes_.size());
    transferSum += probabilities[e] * _PT[e];
  }
  transferSum /= speciesNodes_.size();
}

ScaledValue UndatedDTLModel::getCorrectedTransferExtinctionSum(unsigned int speciesNode) const
{
  return _transferExtinctionSum - _ancestralExctinctionCorrection[speciesNode];
}

ScaledValue UndatedDTLModel::getCorrectedTransferSum(unsigned int geneId, unsigned int speciesId) const
{
  return _survivingTransferSums[geneId] - _ancestralCorrection[geneId][speciesId];
}

void UndatedDTLModel::setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates)
{
  geneRoot_ = 0;
  assert(speciesNodesCount_ == dupRates.size());
  assert(speciesNodesCount_ == lossRates.size());
  assert(speciesNodesCount_ == transferRates.size());
  _PD = dupRates;
  _PL = lossRates;
  _PT = transferRates;
  _PS = std::vector<double>(speciesNodesCount_, 1.0);
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    auto sum = _PD[e] + _PL[e] + _PT[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] /= sum;
  } 
  _uE = std::vector<ScaledValue>(speciesNodesCount_);
  resetTransferSums(_transferExtinctionSum, _ancestralExctinctionCorrection);
  for (unsigned int it = 0; it < IT; ++it) {
    for (auto speciesNode: speciesNodes_) {
      auto e = speciesNode->node_index;
      ScaledValue proba(_PL[e]);
      proba += _uE[e] * _uE[e] * _PD[e] + getCorrectedTransferExtinctionSum(e) * _uE[e];
      if (speciesNode->left) {
        proba += _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index] * _PS[e];
      }
      ASSERT_PROBA(proba)
      _uE[speciesNode->node_index] = proba;
    }
    updateTransferSums(_transferExtinctionSum, _ancestralExctinctionCorrection, _uE);
  }
  invalidateAllCLVs();
}

UndatedDTLModel::~UndatedDTLModel() { }


void UndatedDTLModel::resetTransferSums(ScaledValue &transferSum,
    std::vector<ScaledValue> &ancestralCorrection)
{
  transferSum = ScaledValue();
  ancestralCorrection = std::vector<ScaledValue>(speciesNodesCount_);
}


void UndatedDTLModel::updateCLV(pll_unode_t *geneNode)
{
  auto gid = geneNode->node_index;
  for (auto speciesNode: speciesNodes_) {
    _uq[gid][speciesNode->node_index] = ScaledValue();
  }
  resetTransferSums(_survivingTransferSums[gid], _ancestralCorrection[gid]);
  for (unsigned int it = 0; it < IT; ++it) {
    for (auto speciesNode: speciesNodes_) {
      computeProbability(geneNode, 
          speciesNode, 
          _uq[gid][speciesNode->node_index]);
    }
    updateTransferSums(_survivingTransferSums[gid], _ancestralCorrection[gid], _uq[gid]);
  }
}

pll_rnode_t *UndatedDTLModel::getBestTransfer(unsigned int gid, pll_rnode_t *speciesNode) 
{
  std::unordered_set<unsigned int> parents;
  for (auto parent = speciesNode; parent->parent != 0; parent = parent->parent) {
    parents.insert(parent->node_index);
  }
  pll_rnode_t *bestSpecies = 0;
  ScaledValue bestProba;
  for (auto node: speciesNodes_) {
    if (parents.count(node->node_index)) {
      continue;
    }
    if (bestProba < _uq[gid][node->node_index]) {
      bestProba = _uq[gid][node->node_index];
      bestSpecies = node;
    }
  } 
  return bestSpecies;
}

void UndatedDTLModel::backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot) 
{
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  
  if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
    scenario.addEvent(Scenario::None, gid, e);
    return;
  }
  
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  if (!isGeneLeaf) {
    leftGeneNode = getLeft(geneNode, isVirtualRoot);
    rightGeneNode = getRight(geneNode, isVirtualRoot);
  }
  unsigned int f = 0;
  unsigned int g = 0;
  unsigned int u_left = 0;
  unsigned int u_right = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  
  std::vector<ScaledValue> values(8);
  if (not isGeneLeaf) {
    // S event
    u_left = leftGeneNode->node_index;
    u_right = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      values[0] = _uq[u_left][f] * _uq[u_right][g] * _PS[e];
      values[1] = _uq[u_left][g] * _uq[u_right][f] * _PS[e];
    }
    // D event
    values[2] = _uq[u_left][e];
    values[2] *= _uq[u_right][e];
    values[2] *= _PD[e];
    // T event
    // @todobenoit we should actually look at the max value and not at the sum here
    values[3] = getCorrectedTransferSum(u_left, e) * _uq[u_right][e]; 
    values[4] = getCorrectedTransferSum(u_right, e) * _uq[u_left][e]; 
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[5] = _uq[gid][f] * _uE[g] * _PS[e];
    values[6] = _uq[gid][g] * _uE[f] * _PS[e];
  }
  values[7] = getCorrectedTransferSum(gid, e) * _uE[e];

  unsigned int maxValueIndex = static_cast<unsigned int>(distance(values.begin(), max_element(values.begin(), values.end())));
  // safety check
  if (values[maxValueIndex].isNull()) {
    ScaledValue proba;
    computeProbability(geneNode, speciesNode, proba, isVirtualRoot);
    std::cerr << "warning: null ll scenario " << _uq[gid][e] << " " << proba  << std::endl;
    assert(false);
  }
  pll_rnode_t *dest = 0;
  switch(maxValueIndex) {
    case 0: 
      scenario.addEvent(Scenario::S, gid, e);
      backtrace(leftGeneNode, speciesNode->left, scenario); 
      backtrace(rightGeneNode, speciesNode->right, scenario); 
      break;
    case 1:
      scenario.addEvent(Scenario::S, gid, e);
      backtrace(leftGeneNode, speciesNode->right, scenario); 
      backtrace(rightGeneNode, speciesNode->left, scenario); 
      break;
    case 2:
      scenario.addEvent(Scenario::D, gid, e);
      backtrace(leftGeneNode, speciesNode, scenario); 
      backtrace(rightGeneNode, speciesNode, scenario); 
      break;
    case 3:
      dest =  getBestTransfer(u_left, speciesNode);
      scenario.addTransfer(Scenario::T, gid, e, e, dest->node_index);
      backtrace(leftGeneNode, dest, scenario);
      backtrace(rightGeneNode, speciesNode, scenario);
      break;
    case 4:
      dest =  getBestTransfer(u_right, speciesNode);
      scenario.addTransfer(Scenario::T, gid, e, e, dest->node_index);
      backtrace(rightGeneNode, dest, scenario);
      backtrace(leftGeneNode, speciesNode, scenario);
      break;
    case 5: 
      scenario.addEvent(Scenario::SL, gid, e);
      backtrace(geneNode, speciesNode->left, scenario); 
      break;
    case 6:
      scenario.addEvent(Scenario::SL, gid, e);
      backtrace(geneNode, speciesNode->right, scenario); 
      break;
    case 7:
      dest = getBestTransfer(gid, speciesNode);
      scenario.addTransfer(Scenario::TL, gid, e, e, dest->node_index); 
      backtrace(geneNode, dest, scenario); 
      break;
    default:
      std::cerr << "event " << maxValueIndex << std::endl;
      Logger::error << "Invalid event in UndatedDLModel::backtrace" << std::endl;
      assert(false);
      break;
  }

}


void UndatedDTLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot)
{
  auto gid = geneNode->node_index;
  auto e = speciesNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  
  if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
    proba = ScaledValue(_PS[e], 0);
    return;
  }
  
  ScaledValue oldProba = proba;
  proba = ScaledValue();
  
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  if (!isGeneLeaf) {
    leftGeneNode = getLeft(geneNode, isVirtualRoot);
    rightGeneNode = getRight(geneNode, isVirtualRoot);
  }
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  if (not isGeneLeaf) {
    // S event
    auto u_left = leftGeneNode->node_index;
    auto u_right = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      proba += ScaledValue::superMult1(_uq[u_left][f], _uq[u_right][g],
          _uq[u_left][g], _uq[u_right][f],
          _PS[e]);
    }
    // D event
    ScaledValue temp = _uq[u_left][e];
    temp *= _uq[u_right][e];
    temp *= _PD[e];
    proba += temp;
    // T event
    proba += getCorrectedTransferSum(u_left, e) * _uq[u_right][e]; 
    proba += getCorrectedTransferSum(u_right, e) * _uq[u_left][e]; 
  }
  if (not isSpeciesLeaf) {
    // SL event
    proba += ScaledValue::superMult2(
        _uq[gid][f], _uE[g],
        _uq[gid][g], _uE[f],
        _PS[e]);
  }
  // TL event
  proba += oldProba * getCorrectedTransferExtinctionSum(e);
  proba += getCorrectedTransferSum(gid, e) * _uE[e];

  // DL event
  proba += oldProba * _uE[e] * (2.0 * _PD[e]); 
  //assert(proba.isProba());
}


void UndatedDTLModel::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;;
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    _uq[u][e] = ScaledValue();
  }
  resetTransferSums(_survivingTransferSums[u], _ancestralCorrection[u]);
  for (unsigned int it = 0; it < IT; ++it) {
    for (auto speciesNode: speciesNodes_) {
      unsigned int e = speciesNode->node_index;
      computeProbability(virtualRoot, speciesNode, _uq[u][e], true);
    }
    updateTransferSums(_survivingTransferSums[u], _ancestralCorrection[u], _uq[u]);
  }
}


ScaledValue UndatedDTLModel::getRootLikelihood(pll_unode_t *root) const
{
  ScaledValue sum;
  auto u = root->node_index + _maxGeneId + 1;;
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    sum += _uq[u][e];
  }
  return sum;
}

