#include "UndatedDTLModel.hpp"
#include <Arguments.hpp>
#include <Logger.hpp>

using namespace std;
const int CACHE_SIZE = 100000;
const int IT = 4;

UndatedDTLModel::UndatedDTLModel()
{
  _maxGeneId = 1;
  Logger::info << "creating undated dl model" << endl;
}

#define IS_PROBA(x) ((x) >= 0 && (x) <= 1 && !isnan(x))
#define ASSERT_PROBA(x) assert(IS_PROBA(x));

void UndatedDTLModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  AbstractReconciliationModel::setInitialGeneTree(treeinfo);
  getIdsPostOrder(*treeinfo, _geneIds); 
  _maxGeneId = 0;
  for (auto gid: _geneIds)
    _maxGeneId = max(_maxGeneId, gid);
  // init ua with zeros
  vector<ScaledValue> zeros(speciesNodesCount_);
  _uq = vector<vector<ScaledValue> >(2 * (_maxGeneId + 1),zeros);
  _survivingTransferSums = vector<ScaledValue>(2 * (_maxGeneId + 1));
  _ancestralCorrection = vector<vector<ScaledValue> >(2 * (_maxGeneId + 1),zeros);
  invalidateAllCLVs();
}

void UndatedDTLModel::updateTransferSums(int gid)
{
  // sums
  _survivingTransferSums[gid] = ScaledValue();
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    _survivingTransferSums[gid] += _uq[gid][e] * _PT[e];
  }
  _survivingTransferSums[gid] /= double(speciesNodes_.size());
  
  // ancestral correction
  for (int i = speciesNodes_.size() - 1; i >= 0; --i) {
    auto speciesNode = speciesNodes_[i];
    int e = speciesNode->node_index;
    _ancestralCorrection[gid][e] = _uq[gid][e] * _PT[e];
    if (speciesNode->parent) {
      int p = speciesNode->parent->node_index;
      _ancestralCorrection[gid][e] += _ancestralCorrection[gid][p];
    }
  }
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    _ancestralCorrection[gid][e] /= double(speciesNodes_.size());
  }
}

void UndatedDTLModel::updateExtinctionTransferSums(vector<double> &ancestralExtinctionCorrection)
{
  _transferExtinctionSum = 0.0;
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    _transferExtinctionSum += _PT[e] * _uE[e];
  }
  _transferExtinctionSum /= speciesNodes_.size();

  for (int i = speciesNodes_.size() - 1; i >= 0; --i) {
    auto speciesNode = speciesNodes_[i];
    int e = speciesNode->node_index;
    ancestralExtinctionCorrection[e] = _uE[e] * _PT[e];
    if (speciesNode->parent) {
      int p = speciesNode->parent->node_index;
      ancestralExtinctionCorrection[e] += ancestralExtinctionCorrection[p];
    }
  }
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    ancestralExtinctionCorrection[e] /= double(speciesNodes_.size());
  }
}

void UndatedDTLModel::setRates(double dupRate, 
  double lossRate,
  double transferRate) {
  geneRoot_ = 0;
  _PD = vector<double>(speciesNodesCount_, dupRate);
  _PL = vector<double>(speciesNodesCount_, lossRate);
  _PT = vector<double>(speciesNodesCount_, transferRate);
  _PS = vector<double>(speciesNodesCount_, 1.0);
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    double sum = _PD[e] + _PL[e] + _PT[e] + _PS[e];
    _PD[e] /= sum;
    _PL[e] /= sum;
    _PT[e] /= sum;
    _PS[e] /= sum;
  } 
  _uE = vector<double>(speciesNodesCount_, 0.0);
  vector<double> ancestralExctinctionCorrection(speciesNodesCount_, 0.0);
  _transferExtinctionSum = 0.0;
  for (int it = 0; it < IT; ++it) {
    for (auto speciesNode: speciesNodes_) {
      int e = speciesNode->node_index;
      double proba = _PL[e] + _PD[e] * _uE[e] * _uE[e] + (_transferExtinctionSum - ancestralExctinctionCorrection[e]);
      if (speciesNode->left) {
        proba += _PS[e] * _uE[speciesNode->left->node_index]  * _uE[speciesNode->right->node_index];
      }
      ASSERT_PROBA(proba)
      _uE[speciesNode->node_index] = proba;
    }
    updateExtinctionTransferSums(ancestralExctinctionCorrection);
  }
  invalidateAllCLVs();
}

UndatedDTLModel::~UndatedDTLModel() { }

void UndatedDTLModel::markInvalidatedNodesRec(pll_unode_t *node)
{
  _isCLVUpdated[node->node_index] = false;
  if (node->back->next) {
    markInvalidatedNodesRec(node->back->next);
    markInvalidatedNodesRec(node->back->next->next);
  }
}

void UndatedDTLModel::markInvalidatedNodes(pllmod_treeinfo_t &treeinfo)
{
  for (int nodeIndex: _invalidatedNodes) {
    auto node = treeinfo.subnodes[nodeIndex];
    markInvalidatedNodesRec(node);
  }
  _invalidatedNodes.clear();
}

void UndatedDTLModel::updateCLVsRec(pll_unode_t *node)
{
  if (_isCLVUpdated[node->node_index]) {
    return;
  }
  if (node->next) {
    updateCLVsRec(node->next->back);
    updateCLVsRec(node->next->next->back);
  }
  updateCLV(node);
  _isCLVUpdated[node->node_index] = true;
}

void UndatedDTLModel::updateCLVs(pllmod_treeinfo_t &treeinfo)
{
  markInvalidatedNodes(treeinfo);
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, _geneIds);
  for (auto root: roots) {
    updateCLVsRec(root);
    updateCLVsRec(root->back);
  }
}

void UndatedDTLModel::resetTransferSums(int gid) 
{
  _survivingTransferSums[gid] = ScaledValue();
  for (auto speciesNode: speciesNodes_) {
    int e = speciesNode->node_index;
    _ancestralCorrection[gid][e] = ScaledValue();
  }
}

void UndatedDTLModel::computeGeneProbabilities(pll_unode_t *geneNode)
{
  int gid = geneNode->node_index;
  for (auto speciesNode: speciesNodes_) {
    _uq[gid][speciesNode->node_index] = ScaledValue();
  }
  resetTransferSums(gid);
  for (int it = 0; it < IT; ++it) {
    for (auto speciesNode: speciesNodes_) {
      computeProbability(geneNode, 
          speciesNode, 
          _uq[gid][speciesNode->node_index]);
    }
    updateTransferSums(gid);
  }
}


void UndatedDTLModel::updateCLV(pll_unode_t *geneNode)
{
  computeGeneProbabilities(geneNode);
}

void UndatedDTLModel::invalidateCLV(int nodeIndex)
{
  _invalidatedNodes.insert(nodeIndex);
}
  
void UndatedDTLModel::invalidateAllCLVs()
{
  _isCLVUpdated = vector<bool>(_maxGeneId + 1, false);
}

void UndatedDTLModel::computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot) const
{
  int gid = geneNode->node_index;
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  int e = speciesNode->node_index;
  
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
  int f = 0;
  int g = 0;
  if (!isSpeciesLeaf) {
    f = speciesNode->left->node_index;
    g = speciesNode->right->node_index;
  }
  if (not isGeneLeaf) {
    // S event
    int gp_i = leftGeneNode->node_index;
    int gpp_i = rightGeneNode->node_index;
    if (not isSpeciesLeaf) {
      proba += ScaledValue::superMult1(_uq[gp_i][f], _uq[gpp_i][g],
          _uq[gp_i][g], _uq[gpp_i][f],
          _PS[e]);
    }
    // D event
    ScaledValue temp = _uq[gp_i][e];
    temp *= _uq[gpp_i][e];
    temp *= _PD[e];
    proba += temp;
    // T event
    proba += (_survivingTransferSums[gp_i] - _ancestralCorrection[gp_i][e]) * _uq[gpp_i][e]; 
    proba += (_survivingTransferSums[gpp_i] - _ancestralCorrection[gpp_i][e]) * _uq[gp_i][e]; 
  }
  if (not isSpeciesLeaf) {
    // SL event
    proba += ScaledValue::superMult2(
        _uq[gid][f], _uE[g],
        _uq[gid][g], _uE[f],
        _PS[e]);
  }
  // TL event
  proba += oldProba * _transferExtinctionSum;
  proba += _survivingTransferSums[gid] * _uE[e];

  // DL event
  proba += oldProba * (2.0 * _PD[e] * _uE[e]); 
  //assert(proba.isProba());
}

void UndatedDTLModel::computeLikelihoods(pllmod_treeinfo_t &treeinfo)
{
  vector<ScaledValue> zeros(speciesNodesCount_); 
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, _geneIds);
  for (auto root: roots) {
    int u = root->node_index + _maxGeneId + 1;
    int uprime = root->back->node_index + _maxGeneId + 1;
    for (auto speciesNode: speciesNodes_) {
      int e = speciesNode->node_index;
      _uq[u][e] = ScaledValue();
      _uq[uprime][e] = ScaledValue();
    }
    _survivingTransferSums[u] = ScaledValue();
    for (int it = 0; it < IT; ++it) {
      for (auto speciesNode: speciesNodes_) {
        int e = speciesNode->node_index;
        pll_unode_t virtual_root;
        virtual_root.next = root;
        virtual_root.node_index = root->node_index + _maxGeneId + 1;
        computeProbability(&virtual_root, speciesNode, _uq[u][e], true);
        _uq[uprime][e] = _uq[u][e];
      }
      _survivingTransferSums[u] = ScaledValue();
      for (auto speciesNode: speciesNodes_) {
        int e = speciesNode->node_index;
         _survivingTransferSums[u] += _uq[u][e] * _PT[e];
      }
      _survivingTransferSums[u] /= double(speciesNodes_.size());
      _survivingTransferSums[uprime] = _survivingTransferSums[u];
    }
  }
}

void UndatedDTLModel::updateRoot(pllmod_treeinfo_t &treeinfo) 
{
  vector<pll_unode_t *> roots;
  getRoots(treeinfo, roots, _geneIds);
  // find the best root
  ScaledValue max;
  for (auto root: roots) {
    ScaledValue sum;
    int u = root->node_index + _maxGeneId + 1;;
    for (auto speciesNode: speciesNodes_) {
      int e = speciesNode->node_index;
      sum += _uq[u][e];
    }
    if (max < sum) {
      setRoot(root);
      max = sum;
    }
  }
}
double UndatedDTLModel::getSumLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  ScaledValue total;
  vector<pll_unode_t *> roots;
  getRoots(*treeinfo, roots, _geneIds);
  for (auto root: roots) {
    int u = root->node_index + _maxGeneId + 1;
    ScaledValue sum;
    for (auto speciesNode: speciesNodes_) {
      int e = speciesNode->node_index;
      total += _uq[u][e];
      sum += _uq[u][e];
    }
  }
  double res = total.getLogValue(); //log(root_sum / survive / O_norm * (speciesNodesCount_));
  return res;
}

double UndatedDTLModel::computeLogLikelihoodInternal(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  auto root = getRoot();
  getIdsPostOrder(*treeinfo, _geneIds); 
  updateCLVs(*treeinfo);
  computeLikelihoods(*treeinfo);
  if (Arguments::rootedGeneTree) {// && !getRoot()) {
    updateRoot(*treeinfo);
    while (root != getRoot()) {
      updateCLVs(*treeinfo);
      computeLikelihoods(*treeinfo);
      root = getRoot();
      updateRoot(*treeinfo);
    }
  }
  return getSumLikelihood(treeinfo);
}


