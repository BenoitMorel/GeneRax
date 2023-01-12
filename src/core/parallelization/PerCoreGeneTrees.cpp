#include "PerCoreGeneTrees.hpp"
#include <parallelization/ParallelContext.hpp>
#include <ccp/ConditionalClades.hpp>
#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/LibpllParsers.hpp>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <trees/PLLRootedTree.hpp>

template <typename T>
std::vector<size_t> sort_indexes_descending(const std::vector<T> &v) {
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
  return idx;
}


static std::vector<size_t> getMyIndices(const std::vector<unsigned int> &treeSizes) 
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

static void splitLines(const std::string &input, 
    std::vector<std::string> &output)
{
  std::stringstream ss(input);
  std::string to;
  while(std::getline(ss,to,'\n')){
    if (to.size() > 2) {
      output.push_back(to);
    }
  }
}

std::vector<unsigned int> getCCPSizes(const Families &families) 
{
  unsigned int treesNumber = static_cast<unsigned int>(families.size());
  std::vector<unsigned int> localTreeSizes((treesNumber - 1 ) / ParallelContext::getSize() + 1, 0);
  for (auto i = ParallelContext::getBegin(treesNumber); i < ParallelContext::getEnd(treesNumber); i ++) {
    ConditionalClades cc;
    cc.unserialize(families[i].ccp);
    localTreeSizes[i - ParallelContext::getBegin(treesNumber)] = cc.getCladesNumber();
  }
  std::vector<unsigned int> treeSizes;
  ParallelContext::concatenateUIntVectors(localTreeSizes, treeSizes);
  treeSizes.erase(remove(treeSizes.begin(), treeSizes.end(), 0), treeSizes.end());
  assert(treeSizes.size() == families.size());
  return treeSizes;

}

PerCoreGeneTrees::PerCoreGeneTrees(const Families &families,
    bool acceptMultipleTrees,
    bool ccpMode)
{
  Logger::timed << "Building PerCoreGeneTrees..." << std::endl;
  auto treeSizes = ccpMode ? getCCPSizes(families) : LibpllParsers::parallelGetTreeSizes(families);
  auto myIndices = getMyIndices(treeSizes);

  _geneTrees.resize(myIndices.size());
  unsigned int index = 0;
  for (auto i: myIndices) {
    std::string geneTreeStr;
    FileSystem::getFileContent(families[i].startingGeneTree, geneTreeStr);
    std::vector<std::string> geneTreeStrVector;
    splitLines(geneTreeStr, geneTreeStrVector);
    if (acceptMultipleTrees) {
      if (index == 0) {
        _geneTrees.resize(myIndices.size() * geneTreeStrVector.size());
      }
    } else {
      geneTreeStrVector.resize(1);
    }
    for (auto &currentGeneTreeStr: geneTreeStrVector) {
      _geneTrees[index].name = families[i].name;
      _geneTrees[index].familyIndex = i; 
      _geneTrees[index].mapping.fill(families[i].mappingFile, currentGeneTreeStr);
      _geneTrees[index].geneTree = new PLLUnrootedTree(currentGeneTreeStr, false);
      _geneTrees[index].ownTree = true;
      index++;
    }
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
  corax_rtree_destroy(speciesTree, 0);
  ParallelContext::parallelAnd(ok);
  return ok;
}

void PerCoreGeneTrees::getPerCoreFamilies(const Families &allFamilies,
      Families &perCoreFamilies)
{
  perCoreFamilies.clear();
  PerCoreGeneTrees geneTrees(allFamilies);
  for (const auto &geneTree: geneTrees.getTrees()) {
    perCoreFamilies.push_back(allFamilies[geneTree.familyIndex]);
  }
}

