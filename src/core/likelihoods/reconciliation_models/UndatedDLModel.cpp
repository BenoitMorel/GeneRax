#include "UndatedDLModel.hpp"
#include <IO/Logger.hpp>
#include <algorithm>
#include <util/Scenario.hpp>

#include <cmath>

UndatedDLModel::UndatedDLModel()
{
  _maxGeneId = 1;
}

#define IS_PROBA(x) ((x) >= 0 && (x) <= 1 && !std::isnan(x))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));

void UndatedDLModel::setInitialGeneTree(pll_utree_t *tree)
{
  AbstractReconciliationModel::setInitialGeneTree(tree);
  std::vector<ScaledValue> zeros(speciesNodesCount_);
  _uq = std::vector<std::vector<ScaledValue> >(2 * (_maxGeneId + 1),zeros);
}

static double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

void UndatedDLModel::setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates)
{
  assert(speciesNodesCount_ == dupRates.size());
  assert(speciesNodesCount_ == lossRates.size());
  (void)transferRates;  
  _PD = dupRates;
  _PL = lossRates;
  _PS = std::vector<double>(speciesNodesCount_, 1.0);
  geneRoot_ = 0;
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  _uE = std::vector<double>(speciesNodesCount_, 0.0);
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    double a = _PD[e];
    double b = -1.0;
    double c = _PL[e];
    if (speciesNode->left) {
      c += _PS[e] * _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index];
    }
    double proba = solveSecondDegreePolynome(a, b, c);
    ASSERT_PROBA(proba)
    _uE[speciesNode->node_index] = proba;
  }
  invalidateAllCLVs();
}

UndatedDLModel::~UndatedDLModel() { }

void UndatedDLModel::updateCLV(pll_unode_t *geneNode)
{
  assert(geneNode);
  for (auto speciesNode: speciesNodes_) {
    computeProbability(geneNode, 
        speciesNode, 
        _uq[geneNode->node_index][speciesNode->node_index]);
  }
}


void UndatedDLModel::backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot) 
{
  assert(geneNode);
  assert(speciesNode);
  auto gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;   
  std::vector<ScaledValue> values(5);
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = getLeft(geneNode, isVirtualRoot);
    rightGeneNode = getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !speciesNode->left;
  auto e = speciesNode->node_index;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
    // present
    scenario.addEvent(Scenario::None, gid, e);
    return;
  }
  if (not isGeneLeaf) {
    auto gp_i = leftGeneNode->node_index;
    auto gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      values[0] = _uq[gp_i][f] * _uq[gpp_i][g] * _PS[e];
      values[1] = _uq[gp_i][g] * _uq[gpp_i][f] * _PS[e];
    }
    // D event
    values[2] = _uq[gp_i][e];
    values[2] *= _uq[gpp_i][e];
    values[2] *= _PD[e];
  }
  if (not isSpeciesLeaf) {
    // SL event
    values[3] = _uq[gid][f] * (_uE[g] * _PS[e]);
    values[4] = _uq[gid][g] * (_uE[f] * _PS[e]);
  }

  unsigned maxValueIndex = static_cast<unsigned int>(distance(values.begin(), max_element(values.begin(), values.end())));
  // safety check
  if (values[maxValueIndex] == 0.0) {
    ScaledValue proba = 0.0;
    computeProbability(geneNode, speciesNode, proba, isVirtualRoot);
    std::cerr << "warning: null ll scenario " << _uq[gid][e] << " " << proba  << std::endl;
    assert(false);
  }
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
      scenario.addEvent(Scenario::SL, gid, e);
      backtrace(geneNode, speciesNode->left, scenario); 
      break;
    case 4:
      scenario.addEvent(Scenario::SL, gid, e);
      backtrace(geneNode, speciesNode->right, scenario); 
      break;
    default:
      std::cerr << "event " << maxValueIndex << std::endl;
      Logger::error << "Invalid event in UndatedDLModel::backtrace" << std::endl;
      assert(false);
      break;
  }
}

void UndatedDLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot) const
{
  auto gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = getLeft(geneNode, isVirtualRoot);
    rightGeneNode = getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !speciesNode->left;
  auto e = speciesNode->node_index;
  unsigned int f = 0;
  unsigned int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  proba = ScaledValue();
  if (isSpeciesLeaf and isGeneLeaf) {
    // present
    if (e == geneToSpecies_[gid]) {
      proba = _PS[e];
    } else {
      proba = 0.0;
    }
    return;
  }
  if (not isGeneLeaf) {
    auto gp_i = leftGeneNode->node_index;
    auto gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      /*
      proba += ScaledValue::superMult1(_uq[gp_i][f], _uq[gpp_i][g],
          _uq[gp_i][g], _uq[gpp_i][f],
          _PS[e]);
          */
      proba += (_uq[gp_i][f] * _uq[gpp_i][g] + _uq[gp_i][g] * _uq[gpp_i][f]) * _PS[e];
    }
    // D event
    ScaledValue temp = _uq[gp_i][e];
    temp *= _uq[gpp_i][e];
    temp *= _PD[e];
    proba += temp;
  }
  if (not isSpeciesLeaf) {
    // SL event
    /*
    proba += ScaledValue::superMult2(
        _uq[gid][f], _uE[g],
        _uq[gid][g], _uE[f],
        _PS[e]);
    */
    proba += (_uq[gid][f] * _uE[g] + _uq[gid][g] * _uE[f]) * _PS[e];
  }
  // DL event
  proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
  ASSERT_PROBA(proba);
}
  
ScaledValue UndatedDLModel::getRootLikelihood(pll_unode_t *root) const
{
  ScaledValue sum;
  auto u = root->node_index + _maxGeneId + 1;;
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    sum += _uq[u][e];
  }
  return sum;
}

void UndatedDLModel::accountForSpeciesRoot(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  std::vector<ScaledValue> save_uq(_uq[u]);
  auto leftGeneNode = getLeft(virtualRoot, true);
  auto rightGeneNode = getRight(virtualRoot, true);
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    ScaledValue proba;
    // D
      auto gp_i = leftGeneNode->node_index;
      auto gpp_i = rightGeneNode->node_index;
      ScaledValue temp = _uq[gp_i][e];
      temp *= _uq[gpp_i][e];
      temp *= _PD[e];
      proba += temp;
    // no event
    proba += save_uq[e] * (1.0 - _PD[e]);
    // DL
    proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
    
    _uq[u][e] = proba;
  }
}

void UndatedDLModel::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    computeProbability(virtualRoot, speciesNode, _uq[u][e], true);
  }
//  accountForSpeciesRoot(virtualRoot);
}

