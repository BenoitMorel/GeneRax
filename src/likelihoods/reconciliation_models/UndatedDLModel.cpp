#include "UndatedDLModel.hpp"
#include <IO/Arguments.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <Scenario.hpp>

using namespace std;
const int CACHE_SIZE = 100000;

UndatedDLModel::UndatedDLModel()
{
  _maxGeneId = 1;
  Logger::info << "creating undated dl model" << endl;
}

#define IS_PROBA(x) ((x) >= 0 && (x) <= 1 && !isnan(x))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));

void UndatedDLModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  AbstractReconciliationModel::setInitialGeneTree(treeinfo);
  vector<ScaledValue> zeros(speciesNodesCount_);
  _uq = vector<vector<ScaledValue> >(2 * (_maxGeneId + 1),zeros);
}

static double solveSecondDegreePolynome(double a, double b, double c) 
{
  return 2 * c / (-b + sqrt(b * b - 4 * a * c));
}

void UndatedDLModel::setRates(double dupRate, 
  double lossRate,
  double transferRates) {
  geneRoot_ = 0;
  _PD = vector<double>(speciesNodesCount_, dupRate);
  _PL = vector<double>(speciesNodesCount_, lossRate);
  _PS = vector<double>(speciesNodesCount_, 1.0);
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    double sum = _PD[e] + _PL[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PS[e] /= sum;
  }
  _uE = vector<double>(speciesNodesCount_, 0.0);
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
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
  int gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;   
  vector<ScaledValue> values(5);
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = getLeft(geneNode, isVirtualRoot);
    rightGeneNode = getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !speciesNode->left;
  int e = speciesNode->node_index;
  int f = 0;
  int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
    // present
    return;
  }
  if (not isGeneLeaf) {
    int gp_i = leftGeneNode->node_index;
    int gpp_i = rightGeneNode->node_index;
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

  int maxValueIndex = distance(values.begin(), max_element(values.begin(), values.end()));
  if (values[maxValueIndex].isNull()) {
    ScaledValue proba;
    computeProbability(geneNode, speciesNode, proba, isVirtualRoot);

    cerr << "warning: null ll scenario " << _uq[gid][e] << " " << proba  << endl;
    assert(false);
    //return;
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
      cerr << "event " << maxValueIndex << endl;
      Logger::error << "Invalid event in UndatedDLModel::backtrace" << endl;
      assert(false);
      break;
  }
}

void UndatedDLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot) const
{
  int gid = geneNode->node_index;
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = getLeft(geneNode, isVirtualRoot);
    rightGeneNode = getRight(geneNode, isVirtualRoot);
  }
  bool isSpeciesLeaf = !speciesNode->left;
  int e = speciesNode->node_index;
  int f = 0;
  int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  proba = ScaledValue();
  if (isSpeciesLeaf and isGeneLeaf and e == geneToSpecies_[gid]) {
    // present
    proba = ScaledValue(_PS[e], 0);
    return;
  }
  if (not isGeneLeaf) {
    int gp_i = leftGeneNode->node_index;
    int gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      // S event
      proba += ScaledValue::superMult1(_uq[gp_i][f], _uq[gpp_i][g],
          _uq[gp_i][g], _uq[gpp_i][f],
          _PS[e]);
    }
    // D event
    ScaledValue temp = _uq[gp_i][e];
    temp *= _uq[gpp_i][e];
    temp *= _PD[e];
    proba += temp;
  }
  if (not isSpeciesLeaf) {
    // SL event
    proba += ScaledValue::superMult2(
        _uq[gid][f], _uE[g],
        _uq[gid][g], _uE[f],
        _PS[e]);
  }
  // DL event
  proba /= (1.0 - 2.0 * _PD[e] * _uE[e]); 
  assert(proba.isProba());
}
  
ScaledValue UndatedDLModel::getRootLikelihood(pllmod_treeinfo_t &treeinfo,
    pll_unode_t *root) const
{
  ScaledValue sum;
  int u = root->node_index + _maxGeneId + 1;;
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    sum += _uq[u][e];
  }
  return sum;
}
void UndatedDLModel::computeRootLikelihood(pllmod_treeinfo_t &treeinfo,
    pll_unode_t *virtualRoot)
{
  int u = virtualRoot->node_index;
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    computeProbability(virtualRoot, speciesNode, _uq[u][e], true);
  }
}

