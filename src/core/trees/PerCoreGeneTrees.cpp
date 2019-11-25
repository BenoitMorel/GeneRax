#include "PerCoreGeneTrees.hpp"
#include <parallelization/ParallelContext.hpp>
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


static std::vector<size_t> getMyIndices(const std::vector<unsigned int> &treeSizes, double &loadBalanceRatio) 
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
  double worstLoad = static_cast<double>(*std::max_element(perRankLoad.begin(), perRankLoad.end()));
  loadBalanceRatio = static_cast<double>(averageLoad) / worstLoad;
  return myIndices;
}

PerCoreGeneTrees::PerCoreGeneTrees(const Families &families):
  _loadBalanceRatio(0.0)
{
  auto treeSizes = LibpllParsers::parallelGetTreeSizes(families);
  auto myIndices = getMyIndices(treeSizes, _loadBalanceRatio);

  _geneTrees.resize(myIndices.size());
  unsigned int index = 0;
  for (auto i: myIndices) {
    std::string geneTreeStr;
    FileSystem::getFileContent(families[i].startingGeneTree, geneTreeStr);
    _geneTrees[index].name = families[i].name;
    _geneTrees[index].mapping.fill(families[i].mappingFile, geneTreeStr);
    _geneTrees[index].geneTree = new PLLUnrootedTree(families[i].startingGeneTree, true);
    _geneTrees[index].ownTree = true;
    index++;
  }
  ParallelContext::barrier();
}
  
PerCoreGeneTrees::PerCoreGeneTrees(const GeneSpeciesMapping &mapping, PLLUnrootedTree &geneTree)
{
  if (ParallelContext::getRank() == 0) {
    _geneTrees.resize(1);
    _geneTrees[0].name = "JointTree";
    _geneTrees[0].mapping = mapping;
    _geneTrees[0].geneTree = &geneTree;
    _geneTrees[0].ownTree = false;
  }
  ParallelContext::barrier();
}
  

bool PerCoreGeneTrees::checkMappings(const std::string &speciesTreeFile)
{
  bool ok = true;
  auto *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile);
  for (auto &tree: _geneTrees) {
    if (!tree.mapping.check(tree.geneTree->getRawPtr(), speciesTree)) {
      ok = false;
      Logger::error << "Invalid mapping for tree " << tree.name << std::endl;
    }
  }
  pll_rtree_destroy(speciesTree, 0);
  ParallelContext::parallelAnd(ok);
  return ok;
}


