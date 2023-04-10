#include "TrimFamilies.hpp"

#include <vector>
#include <algorithm>

#include <ccp/ConditionalClades.hpp>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <parallelization/ParallelContext.hpp>
#include <util/types.hpp>

void TrimFamilies::trimHighCladesNumber(Families &families,
    double keepRatio)
{
  std::vector<unsigned int> localCcpSizes;
  std::vector<unsigned int> ccpSizes;
  std::vector<unsigned int> localCheck;
  std::vector<unsigned int>check;
  auto begin = ParallelContext::getBegin(families.size());
  auto end = ParallelContext::getEnd(families.size());
  for (auto i = begin; i < end; ++i) {
    ConditionalClades ccp;
    ccp.unserialize(families[i].ccp);
    localCcpSizes.push_back(ccp.getCladesNumber());
    localCheck.push_back(i);
  }
  ParallelContext::concatenateHetherogeneousUIntVectors(
      localCcpSizes, ccpSizes);
  ParallelContext::concatenateHetherogeneousUIntVectors(
      localCheck, check);
  std::vector<PairUInt> sizeToIndex(families.size());
  for (unsigned int i = 0; i < families.size(); ++i) {
    sizeToIndex[i].first = ccpSizes[i];
    sizeToIndex[i].second = i;
  }
  std::sort(sizeToIndex.begin(), sizeToIndex.end());
  Families copy = families;
  unsigned int cutAfter = (unsigned int)(double(families.size()) * keepRatio);
  families.clear();
  assert(cutAfter > 0);
  for (unsigned int i = 0; i < cutAfter; ++i) {
    assert(ParallelContext::isIntEqual(sizeToIndex[i].second));
    families.push_back(copy[sizeToIndex[i].second]);
  }
  Logger::info << "Trimming families with from " << sizeToIndex[cutAfter].first << " to " << sizeToIndex.back().first << " clades" << std::endl;
}

void TrimFamilies::trimMinSpeciesCoverage(Families &families,
      unsigned int minCoverage)
{
  auto begin = ParallelContext::getBegin(families.size());
  auto end = ParallelContext::getEnd(families.size());
  std::vector<unsigned int> localToKeep;
  std::vector<unsigned int> toKeep;
  for (unsigned int i = begin; i < end; ++i) {
    auto &family = families[i];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    if (mapping.getCoveredSpecies().size() >= minCoverage) {
      localToKeep.push_back(i); 
    }
  }
  ParallelContext::concatenateHetherogeneousUIntVectors(localToKeep,
      toKeep);
  Families copy = families;
  families.clear();
  for (auto i: toKeep) {
    assert(ParallelContext::isIntEqual(i));
    families.push_back(copy[i]);
  }
}
  
void TrimFamilies::trimCladeSplitRatio(Families &families,
      double maxRatio)
{
  auto begin = ParallelContext::getBegin(families.size());
  auto end = ParallelContext::getEnd(families.size());
  std::vector<unsigned int> localToKeep;
  std::vector<unsigned int> toKeep;
  unsigned int allClades = 0;
  unsigned int filteredClades = 0;
  for (unsigned int i = begin; i < end; ++i) {
    ConditionalClades ccp;
    ccp.unserialize(families[i].ccp);
    auto c = ccp.getCladesNumber();
    auto n = ccp.getLeafNumber() * 2;
    allClades += c;
    if (maxRatio * n >= c) {
      localToKeep.push_back(i); 
      filteredClades += c;
    }
  }
  ParallelContext::concatenateHetherogeneousUIntVectors(localToKeep,
      toKeep);
  ParallelContext::sumUInt(allClades);
  ParallelContext::sumUInt(filteredClades);
  Logger::info << "Clades before trimming: " << allClades << std::endl;
  Logger::info << "Clades after trimming: " << filteredClades << std::endl;
  Families copy = families;
  families.clear();
  for (auto i: toKeep) {
    assert(ParallelContext::isIntEqual(i));
    families.push_back(copy[i]);
  }
}




