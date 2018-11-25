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
    double subdivisionSize = 0.05;
    int minSubdivisions = 5;
    if (minSubdivisions * subdivisionSize > branchLength) {
      subdivisionSize = branchLength / minSubdivisions;
    }
    branchSubdivisions[speciesId].push_back(0);
    while (branchLength > subdivisionSize) {
      branchSubdivisions[speciesId].push_back(subdivisionSize);
      branchLength -= subdivisionSize;
    }
    if (branchLength > 0)
      branchSubdivisions[speciesId].push_back(branchLength);
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
    /*
    for (int s = 0; s < subdivisions; ++s) {
      Logger::info << extinctionProba_[speciesId][s] << " ";
    }
    Logger::info << endl;
    */
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
      propagationProba_[speciesId][s] = propagatePropagationProba(extinctionProba_[speciesId][s-1],
          branchSubdivisions_[speciesId][s]);
    }
    /*
    for (int s = 0; s < subdivisions; ++s) {
      Logger::info << propagationProba_[speciesId][s] << " ";  
    }
    Logger::info << endl;
    */
  }
}

double DatedDLModel::propagatePropagationProba(double initialProba, double branchLength)
{
  if (dupRate_ == lossRate_) {
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
  vector<pll_unode_t *> roots;
  getRoots(*treeinfo, roots, geneIds);
  double ll = 0;
  double norm = 0;
  bool selectMax = true;
  for (auto geneRoot: roots) {
    for (auto species: speciesNodes_) {
      int speciesId = species->node_index;
      for (int i = 0; i < branchSubdivisions_[speciesId].size(); ++i) {
        double localLL =  getRecProba(geneRoot->node_index, speciesId, i) * getRecProba(geneRoot->back->node_index, speciesId, i); 
        if (selectMax) {
          ll = max(ll, localLL);
        } else {
          ll += localLL;
          norm +=  1 - extinctionProba_[speciesId][i];
        }
      }
    }
  }
  if (!selectMax) 
    ll /= norm;
  return ll;
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
      clv[speciesId][subdivision] = computeRecProbaIntraBranch(geneNode, speciesNode, subdivision);
    }
    for (int i = 0; i < subdivisions; ++i) {
      double d = clv[speciesId][i];
      if (d < 0 || d > 1.0)
        Logger::info << d << endl;
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
  
double DatedDLModel::computeRecProbaIntraBranch(pll_unode_t *geneNode, pll_rnode_t *speciesNode, int subdivision)
{
  bool isGeneLeaf = !geneNode->next;
  int geneId = geneNode->node_index;
  int speciesId = speciesNode->node_index;
  double res = 0.0;
  // No event case
  res += propagationProba_[speciesId][subdivision] * getRecProba(geneId, speciesId, subdivision - 1);
  // duplication case
  if (!isGeneLeaf) {
    int leftGeneId = geneNode->next->back->node_index;
    int rightGeneId = geneNode->next->next->back->node_index;
    double l = branchSubdivisions_[speciesId][subdivision];
    double leftProba = getRecProba(leftGeneId, speciesId, subdivision - 1);
    double rightProba = getRecProba(rightGeneId, speciesId, subdivision - 1);
    res += 2.0 * dupRate_ * l * leftProba * rightProba;
  }
  return res;

}

double DatedDLModel::getRecProba(int geneId, int speciesId)
{
  return clvs_[geneId].clv[speciesId].back();
}

double DatedDLModel::getRecProba(int geneId, int speciesId, int subdivision)
{
  return clvs_[geneId].clv[speciesId][subdivision];
}

double DatedDLModel::getExtProba(int speciesId)
{
  return extinctionProba_[speciesId].back();
}

