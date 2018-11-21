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
  int speciesNumber = speciesTree->tip_count + speciesTree->inner_count;
  extinctionProba_.resize(speciesNumber);
  for (int i = 0; i < speciesNumber; ++i) {
    auto node = speciesTree->nodes[i];
    int speciesId = node->node_index;
    int subdivisions = branchSubdivisions_[speciesId].size();
    extinctionProba_[speciesId].resize(subdivisions);
    // first branch subdivision
    if (!node->left) {
      // species leaf
      extinctionProba_[speciesId][0] = 1 - probaGeneSampled_; 
    } else {
      // species internal node
      extinctionProba_[speciesId][0] = extinctionProba_[node->left->node_index].back() *
        extinctionProba_[node->right->node_index].back();
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
  int speciesNumber = speciesTree->tip_count + speciesTree->inner_count;
  propagationProba_.resize(speciesNumber);
  for (int i = 0; i < speciesNumber; ++i) {
    auto node = speciesTree->nodes[i];
    int speciesId = node->node_index;
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
  // todobenoit: check that we do not need to check that speciesTree nodes
  // are ordered with postorder traversal (we do it in UndatedDLModel)
  speciesTree_ = speciesTree;
  // build subdivisions
  buildSubdivisions(speciesTree, branchSubdivisions_);
}


void DatedDLModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{

}

void DatedDLModel::setGeneSpeciesMap(const GeneSpeciesMapping &map)
{

}

double DatedDLModel::computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{
  vector<int> geneIds;
  getIdsPostOrder(*treeinfo, geneIds);
  for (auto geneId: geneIds) {
    updateCLV(treeinfo->subnodes[geneIds[geneId]]);
  }
  return 0;
}

void DatedDLModel::updateCLV(pll_unode_t *geneNode)
{
  Logger::info << "update clv " << geneNode->node_index << endl;
}



