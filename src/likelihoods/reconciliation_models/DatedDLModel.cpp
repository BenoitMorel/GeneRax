#include "DatedDLModel.hpp"
#include <Logger.hpp>

static const double EPSILON = 0.0000000001;

void buildSubdivisions(pll_rtree_t *speciesTree,
  vector<vector<double> > &branchSubdivisions) 
{
  int maxSpeciesNodeIndex = speciesTree->tip_count + speciesTree->inner_count;
  branchSubdivisions = vector<vector<double> >(maxSpeciesNodeIndex);
  int subdivisionsNumber = 5;
  for (int i = 0; i < maxSpeciesNodeIndex; ++i) {
    int speciesId = speciesTree->nodes[i]->node_index;
    double branchLength = max(speciesTree->nodes[i]->length, EPSILON);
    // todobenoit: have different number of subdivisions depending on the
    // length of the branchs
    branchSubdivisions[speciesId] = vector<double>(subdivisionsNumber + 1, branchLength / subdivisionsNumber);
    branchSubdivisions[speciesId][0] = 0.0;
  }
}

void DatedDLModel::computeExtinctionProbas(pll_rtree_t *speciesTree)
{
  extinctionProba_.resize(speciesNodesCount_);
  for (auto speciesNode: speciesNodes_) {
    int speciesId = speciesNode->node_index;
    int subdivisions = branchSubdivisions_[speciesId].size();
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
    for (int s = 1; s < subdivisions; ++s) {
      extinctionProba_[speciesId][s] = propagateExtinctionProba(extinctionProba_[speciesId][s-1],
          branchSubdivisions_[speciesId][s]);
    }
    
  }
}

double DatedDLModel::propagateExtinctionProba(double initialProba, double branchLength)
{
  if (dupRate_ == lossRate_) {
    double denom = lossRate_ * (initialProba - 1.0) * branchLength - 1.0;
    return 1 + (1 - initialProba) / denom;
  }
  double a = exp(diffRates_ * branchLength);
  double b = dupRate_ * (initialProba - 1.0);
  double c = lossRate_ - dupRate_ * initialProba;
  return (lossRate_ + diffRates_ / (1.0 + a*b/c)) / dupRate_;
}

void DatedDLModel::computePropagationProbas(pll_rtree_t *speciesTree)
{
  propagationProba_.resize(speciesNodesCount_);
  for (auto speciesNode: speciesNodes_) {
    int speciesId = speciesNode->node_index;
    int subdivisions = branchSubdivisions_[speciesId].size();
    propagationProba_[speciesId].resize(subdivisions);
    // first branch subdivision
    propagationProba_[speciesId][0] = 1.0; 
    // go up in the branch
    for (int s = 1; s < subdivisions; ++s) {
      propagationProba_[speciesId][s] = propagatePropagationProba(propagationProba_[speciesId][s-1],
          branchSubdivisions_[speciesId][s]);
    }
    
  }
}

double DatedDLModel::propagatePropagationProba(double initialProba, double branchLength)
{
  if (dupRate_ == lossRate_) {
    double x = lossRate_ * (initialProba - 1.0) * branchLength - 1.0;
    return 1.0 / pow(x, 2.0);
  }
  double x = exp(diffRates_ * branchLength);
  double a = x * pow(diffRates_, 2.0);
  double b = dupRate_ - x * lossRate_;
  double c = (x - 1.0) * dupRate_ * initialProba;
  return a / pow(b + c, 2.0);
}



DatedDLModel::DatedDLModel():
  probaGeneSampled_(1.0)
{

}

void DatedDLModel::setRates(double dupRate, double lossRate, double transferRate)
{
  if (dupRate == 0)
    dupRate = EPSILON;
  if (lossRate == 0)
    lossRate = EPSILON;
  dupRate_ = dupRate;
  lossRate_ = lossRate;
  diffRates_ = dupRate - lossRate;
  
  computeExtinctionProbas(speciesTree_);
  computePropagationProbas(speciesTree_);
}

void DatedDLModel::setSpeciesTree(pll_rtree_t *speciesTree)
{
  AbstractReconciliationModel::setSpeciesTree(speciesTree);
  // todobenoit: check that we do not need to check that speciesTree nodes
  // are ordered with postorder traversal (we do it in UndatedDLModel)
  // build subdivisions
  buildSubdivisions(speciesTree, branchSubdivisions_);
}


double DatedDLModel::computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  vector<int> geneIds;
  getIdsPostOrder(*treeinfo, geneIds);
  clvs_.resize(geneIds.size());
  for (auto geneId: geneIds) {
    updateCLV(treeinfo->subnodes[geneId]);
  }
  return 0;
}

void DatedDLModel::updateCLV(pll_unode_t *geneNode)
{
  int geneId = geneNode->node_index;
  auto &clv = clvs_[geneId].clv;
  clv = vector<vector<double > >(speciesNodesCount_);
  for (auto speciesNode: speciesNodes_) {
    int speciesId = speciesNode->node_index;
    int subdivisions = branchSubdivisions_[speciesId].size();
    clv[speciesId] = vector<double>(subdivisions, 0.0);
    clv[speciesId][0] = computeRecProbaInterBranch(geneNode, speciesNode);
    for (int subdivision = 1; subdivision < subdivisions; ++subdivision) {
      //clv[speciesId][subdivision] = computeRecProbaIntraBranch();
    }
  }
}

double DatedDLModel::computeRecProbaInterBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode)
{
  bool isGeneLeaf = !geneNode->next;
  bool isSpeciesLeaf = !speciesNode->left;
  int geneId = geneNode->node_index;
  int speciesId = speciesNode->node_index;
  
  if (isGeneLeaf && isSpeciesLeaf) {
    // trivial case
    return (geneToSpecies_[geneId] == speciesId) ? 1 : 0; 
  }
  if (isSpeciesLeaf) {
    return 0.0; // todobenoit I am not sure about this
  }
  int leftSpeciesId = speciesNode->left->node_index;
  int rightSpeciesId = speciesNode->right->node_index;

  double res = 0.0;
  // todobenoit check that...
  // Speciation-Loss case
  res += getRecProba(geneId, leftSpeciesId) * getExtProba(rightSpeciesId);  
  res += getRecProba(geneId, rightSpeciesId) * getExtProba(leftSpeciesId);  
  
  if (!isGeneLeaf) {
    int leftGeneId = geneNode->next->back->node_index;
    int rightGeneId = geneNode->next->next->back->node_index;
    // Speciation case
    res += getRecProba(leftGeneId, leftSpeciesId) * getRecProba(rightGeneId, rightSpeciesId);
    res += getRecProba(rightGeneId, leftSpeciesId) * getRecProba(leftGeneId, rightSpeciesId);

  }
  return res;
}

double DatedDLModel::getRecProba(int geneId, int speciesId)
{
  return clvs_[geneId].clv[speciesId].back();
}

double DatedDLModel::getExtProba(int speciesId)
{
  return extinctionProba_[speciesId].back();
}

