#include "DatedDLModel.hpp"
#include <iostream>


void buildSubdivisions(pll_rtree_t *speciesTree,
  vector<vector<double> > &branchSubdivisions) 
{
  int maxSpeciesNodeIndex = speciesTree->tip_count + speciesTree->inner_count;
  branchSubdivisions = vector<vector<double> >(maxSpeciesNodeIndex);
  for (int i = 0; i < maxSpeciesNodeIndex; ++i) {
    int speciesId = speciesTree->nodes[i]->node_index;
    double branchLength = speciesTree->nodes[i]->length;
    cout << "BL: " << branchLength << endl;
  }
}



void DatedDLModel::setRates(double dupRate, double lossRate, double transferRate)
{

}

void DatedDLModel::setSpeciesTree(pll_rtree_t *speciesTree)
{
  // build subdivisions
  buildSubdivisions(speciesTree, branchSubdivisions_);
  // compute extinctionProba_

  // compute propagationProba_
}


void DatedDLModel::setInitialGeneTree(shared_ptr<pllmod_treeinfo_t> treeinfo)
{

}

void DatedDLModel::setGeneSpeciesMap(const GeneSpeciesMapping &map)
{

}

void DatedDLModel::setRoot(pll_unode_t * root)
{
  geneRoot = root;
}

pll_unode_t *DatedDLModel::getRoot()
{
  return geneRoot;
}

double DatedDLModel::computeLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo)
{

  return 0;
}


