#include "PerCoreGeneTrees.hpp"
#include <ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/LibpllParsers.hpp>
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
    currentRank = (currentRank + 1) % ParallelContext::getSize();
    for (; perRankLoad[currentRank] > averageLoad;(currentRank = (currentRank + 1) % ParallelContext::getSize())) {}
  }
  return myIndices;
}

PerCoreGeneTrees::PerCoreGeneTrees(const vector<FamiliesFileParser::FamilyInfo> &families)
{
  vector<int> treeSizes = LibpllParsers::parallelGetTreeSizes(families);
  vector<size_t> myIndices = getMyIndices(treeSizes);

  for (auto i: myIndices) {
    string geneTreeStr;
    FileSystem::getFileContent(families[i].startingGeneTree, geneTreeStr);
    _geneTrees.push_back(GeneTree());
    _geneTrees.back().mapping.fill(families[i].mappingFile, geneTreeStr);
    _geneTrees.back().tree = LibpllParsers::readNewickFromFile(families[i].startingGeneTree);
  }
}
