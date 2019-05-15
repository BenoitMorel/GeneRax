#include "PerCoreGeneTrees.hpp"
#include <ParallelContext.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <numeric>
#include <iostream>

template <typename T>
std::vector<size_t> sort_indexes_descending(const std::vector<T> &v) {
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
  return idx;
}


std::vector<size_t> getMyIndices(const std::vector<unsigned int> &treeSizes) 
{
  std::vector<size_t> sortedIndices = sort_indexes_descending<unsigned int>(treeSizes);
  std::vector<size_t> myIndices;
  std::vector<unsigned int> perRankLoad(ParallelContext::getSize(), 0);
  unsigned int averageLoad = 0;
  for (auto load: treeSizes) {
    averageLoad += load;
  }
  averageLoad /= ParallelContext::getSize();
  unsigned int currentRank = 0;
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

PerCoreGeneTrees::PerCoreGeneTrees(const std::vector<FamiliesFileParser::FamilyInfo> &families)
{
  auto treeSizes = LibpllParsers::parallelGetTreeSizes(families);
  auto myIndices = getMyIndices(treeSizes);

  _geneTrees.resize(myIndices.size());
  unsigned int index = 0;
  for (auto i: myIndices) {
    std::string geneTreeStr;
    FileSystem::getFileContent(families[i].startingGeneTree, geneTreeStr);
    _geneTrees[index].name = families[i].name;
    _geneTrees[index].mapping.fill(families[i].mappingFile, geneTreeStr);
    _geneTrees[index].tree = LibpllParsers::readNewickFromFile(families[i].startingGeneTree);
    index++;
  }
}
