#include "PerCoreGeneTrees.hpp"
#include <ParallelContext.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/Logger.hpp>
#include <algorithm>
#include <numeric>
#include <iostream>

template <typename T>
vector<size_t> sort_indexes_descending(const vector<T> &v) {
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
  return idx;
}

vector<int> getTreeSizes(const vector<FamiliesFileParser::FamilyInfo> &families) 
{
  int treesNumber = families.size();
  vector<int> localTreeSizes((treesNumber - 1 ) / ParallelContext::getSize() + 1, 0);
  for (int i = ParallelContext::getBegin(treesNumber); i < ParallelContext::getEnd(treesNumber); i ++) {
    pll_utree_t *tree = LibpllEvaluation::readNewickFromFile(families[i].startingGeneTree);
    int taxa = tree->tip_count;
    localTreeSizes[i - ParallelContext::getBegin(treesNumber)] = taxa;
    pll_utree_destroy(tree, 0);
  }
  vector<int> treeSizes;
  ParallelContext::concatenateIntVectors(localTreeSizes, treeSizes);
  treeSizes.erase(std::remove(treeSizes.begin(), treeSizes.end(), 0), treeSizes.end());
  assert(treeSizes.size() == families.size());
  return treeSizes;
}

vector<size_t> getMyIndices(const vector<int> &treeSizes) 
{
  vector<size_t> sortedIndices = sort_indexes_descending<int>(treeSizes);
  vector<size_t> myIndices;
  vector<int> perRankLoad(ParallelContext::getSize(), 0);
  int averageLoad = 0;
  for (int load: treeSizes) {
    averageLoad += load;
  }
  averageLoad /= ParallelContext::getSize();
  int currentRank = 0;
  for (auto index: sortedIndices) {
    perRankLoad[currentRank] += treeSizes[index];
    if (currentRank == ParallelContext::getRank()) {
      myIndices.push_back(index);
    }
    for (; perRankLoad[currentRank] > averageLoad;(currentRank = (currentRank + 1) % ParallelContext::getSize())) {}
  }
  return myIndices;
}

PerCoreGeneTrees::PerCoreGeneTrees(const vector<FamiliesFileParser::FamilyInfo> &families)
{
  vector<int> treeSizes = getTreeSizes(families);
  vector<size_t> myIndices = getMyIndices(treeSizes);

  for (auto i: myIndices) {
    _geneTrees.push_back(GeneTree());
    _geneTrees.back().mapping = GeneSpeciesMapping(families[i].mappingFile);
    _geneTrees.back().tree = LibpllEvaluation::readNewickFromFile(families[i].startingGeneTree);
  }
}
