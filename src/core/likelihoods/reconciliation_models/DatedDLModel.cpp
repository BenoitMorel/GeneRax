#include "DatedDLModel.hpp"
#include <IO/Logger.hpp>
#include <cmath>

static const double EPSILON = 0.0000000001;

void DatedDLModel::setInitialGeneTree(pll_utree_t *tree)
{
  AbstractReconciliationModel::setInitialGeneTree(tree);
  clvs_ = std::vector<DDL_CLV>((_maxGeneId + 1) * 2); 
}

void buildSubdivisions(pll_rtree_t *speciesTree,
  std::vector<std::vector<double> > &branchSubdivisions) 
{
  auto maxSpeciesNodeIndex = speciesTree->tip_count + speciesTree->inner_count;
  branchSubdivisions = std::vector<std::vector<double> >(maxSpeciesNodeIndex);
  for (unsigned int i = 0; i < maxSpeciesNodeIndex; ++i) {
    auto speciesId = speciesTree->nodes[i]->node_index;
    double branchLength = std::max(speciesTree->nodes[i]->length, 0.0001);
    double subdivisionSize = 0.05;
    unsigned int minSubdivisions = 5;
    if (minSubdivisions * subdivisionSize > branchLength) {
      subdivisionSize = branchLength / minSubdivisions;
    }
    branchSubdivisions[speciesId].push_back(0);
    while (branchLength > subdivisionSize) {
      branchSubdivisions[speciesId].push_back(subdivisionSize);
      branchLength -= subdivisionSize;
    }
    if (branchLength > EPSILON)
      branchSubdivisions[speciesId].push_back(branchLength);
  }
}

void DatedDLModel::computeExtinctionProbas()
{
  extinctionProba_.resize(speciesNodesCount_);
  for (auto speciesNode: speciesNodes_) {
    auto speciesId = speciesNode->node_index;
    unsigned int subdivisions = static_cast<unsigned int>(branchSubdivisions_[speciesId].size());
    extinctionProba_[speciesId].resize(subdivisions);
    // first branch subdivision
    if (!speciesNode->left) {
      // species leaf
      extinctionProba_[speciesId][0] = 1 - probaGeneSampled_; 
    } else {
      // species internal node
      extinctionProba_[speciesId][0] = getExtProba(speciesNode->left->node_index) *
        getExtProba(speciesNode->right->node_index);
    }
    // go up in the branch
    for (unsigned int s = 1; s < subdivisions; ++s) {
      extinctionProba_[speciesId][s] = propagateExtinctionProba(extinctionProba_[speciesId][s-1],
          branchSubdivisions_[speciesId][s]);
    }
  }
}

double DatedDLModel::propagateExtinctionProba(double initialProba, double branchLength)
{
  assert(branchLength > EPSILON);
  if (fabs(diffRates_) < EPSILON) {
    double denom = lossRate_ * (initialProba - 1.0) * branchLength - 1.0;
    return 1 + (1 - initialProba) / denom;
  }
  double a = exp(diffRates_ * branchLength);
  double b = dupRate_ * (initialProba - 1.0);
  double c = lossRate_ - dupRate_ * initialProba;
  double res = (lossRate_ + diffRates_ / (1.0 + a*b/c)) / dupRate_;
  assert(std::isnormal(res));
  return res;
}

void DatedDLModel::computePropagationProbas()
{
  propagationProba_.resize(speciesNodesCount_);
  for (auto speciesNode: speciesNodes_) {
    auto speciesId = speciesNode->node_index;
    auto subdivisions = branchSubdivisions_[speciesId].size();
    propagationProba_[speciesId].resize(subdivisions);
    // first branch subdivision
    propagationProba_[speciesId][0] = 1.0; 
    // go up in the branch
    for (unsigned int s = 1; s < subdivisions; ++s) {
      propagationProba_[speciesId][s] = propagatePropagationProba(extinctionProba_[speciesId][s-1],
          branchSubdivisions_[speciesId][s]);
      assert(std::isnormal(propagationProba_[speciesId][s]));
    }
  }
}

double DatedDLModel::propagatePropagationProba(double initialProba, double branchLength)
{
  assert(branchLength > EPSILON);
  if (fabs(diffRates_) < EPSILON) {
    double x = lossRate_ * (initialProba - 1.0) * branchLength - 1.0;
    double res = 1.0 / pow(x, 2.0);
    return res;
  }
  double x = exp(-diffRates_ * branchLength);
  double a = x * pow(diffRates_, 2.0);
  double b = dupRate_ - x * lossRate_;
  double c = (x - 1.0) * dupRate_ * initialProba;
  double res = a / pow(b + c, 2.0);
  return res;
}



DatedDLModel::DatedDLModel():
  probaGeneSampled_(1.0)
{
}

void DatedDLModel::setRates(double dupRate, double lossRate, double transferRate, const std::vector<double> &speciesScalers)
{
  (void)transferRate;
  (void)speciesScalers;
  if (0.0 == dupRate)
    dupRate = EPSILON;
  if (0.0 == lossRate)
    lossRate = EPSILON;
  dupRate_ = dupRate;
  lossRate_ = lossRate;
  diffRates_ = dupRate - lossRate;
  computeExtinctionProbas();
  computePropagationProbas();
  invalidateAllCLVs();
}

void DatedDLModel::setSpeciesTree(pll_rtree_t *speciesTree)
{
  AbstractReconciliationModel::setSpeciesTree(speciesTree);
  buildSubdivisions(speciesTree, branchSubdivisions_);
}


void DatedDLModel::updateCLV(pll_unode_t *geneNode)
{
  unsigned int geneId = geneNode->node_index;
  auto &clv = clvs_[geneId].clv;
  clv = std::vector<std::vector<ScaledValue> >(speciesNodesCount_);
  for (auto speciesNode: speciesNodes_) {
    computeCLVCell(geneNode, speciesNode, clv[speciesNode->node_index], false);
  }
}

void DatedDLModel::computeCLVCell(pll_unode_t *geneNode, pll_rnode_t *speciesNode, std::vector<ScaledValue> &speciesCell, bool isVirtualRoot)
{
  unsigned int subdivisions = static_cast<unsigned int>(branchSubdivisions_[speciesNode->node_index].size());
  speciesCell = std::vector<ScaledValue>(subdivisions);
  speciesCell[0] = computeRecProbaInterBranch(geneNode, speciesNode, isVirtualRoot);
  for (unsigned int subdivision = 1; subdivision < subdivisions; ++subdivision) {
    speciesCell[subdivision] = computeRecProbaIntraBranch(geneNode, speciesNode, subdivision, isVirtualRoot);
  }
}

ScaledValue DatedDLModel::computeRecProbaInterBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode, bool isVirtualRoot)
{
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  auto geneId = geneNode->node_index;
  auto speciesId = speciesNode->node_index;
  
  if (isGeneLeaf && isSpeciesLeaf) {
    // trivial case
    return (geneToSpecies_[geneId] == speciesId) ? ScaledValue(1) : ScaledValue(); 
  }
  ScaledValue res;
  if (isSpeciesLeaf) {
    return res; // todobenoit I am not sure about this
  }
  auto leftSpeciesId = speciesNode->left->node_index;
  auto rightSpeciesId = speciesNode->right->node_index;


  res += getRecProba(geneId, leftSpeciesId) * getExtProba(rightSpeciesId);  
  res += getRecProba(geneId, rightSpeciesId) * getExtProba(leftSpeciesId);  
  if (!isGeneLeaf) {
    auto leftGeneId = getLeft(geneNode, isVirtualRoot)->node_index;
    auto rightGeneId = getRight(geneNode, isVirtualRoot)->node_index;
    // Speciation case
    res += getRecProba(leftGeneId, leftSpeciesId) * getRecProba(rightGeneId, rightSpeciesId);
    res += getRecProba(rightGeneId, leftSpeciesId) * getRecProba(leftGeneId, rightSpeciesId);
  }
  return res;
}
  
ScaledValue DatedDLModel::computeRecProbaIntraBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode, unsigned int subdivision, bool isVirtualRoot)
{
  bool isGeneLeaf = !geneNode->next;
  auto geneId = geneNode->node_index;
  auto speciesId = speciesNode->node_index;
  ScaledValue res;
  // No event case
  res += getRecProba(geneId, speciesId, subdivision - 1) * propagationProba_[speciesId][subdivision];
  // duplication case
  if (!isGeneLeaf) {
    auto leftGeneId = getLeft(geneNode, isVirtualRoot)->node_index;
    auto rightGeneId = getRight(geneNode, isVirtualRoot)->node_index;
    double l = branchSubdivisions_[speciesId][subdivision];
    auto leftProba = getRecProba(leftGeneId, speciesId, subdivision - 1); 
    auto rightProba = getRecProba(rightGeneId, speciesId, subdivision - 1);
    res += leftProba * rightProba * dupRate_ * l;
  }
  return res;

}


double DatedDLModel::getExtProba(unsigned int speciesId)
{
  return extinctionProba_[speciesId].back();
}

ScaledValue DatedDLModel::getRootLikelihood(pll_unode_t *root) const
{
  ScaledValue sum;
  auto u = root->node_index + _maxGeneId + 1;;
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    sum += clvs_[u].clv[e].back();
  }
  return sum;
}

void DatedDLModel::computeRootLikelihood(pll_unode_t *virtualRoot)
{
  auto u = virtualRoot->node_index;
  clvs_[u].clv = std::vector<std::vector<ScaledValue> >(speciesNodesCount_);
  for (auto speciesNode: speciesNodes_) {
    auto e = speciesNode->node_index;
    computeCLVCell(virtualRoot, speciesNode, clvs_[u].clv[e], true);
  }
}


