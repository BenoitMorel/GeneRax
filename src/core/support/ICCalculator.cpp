#include "ICCalculator.hpp"

#include <array>
#include <unordered_map>
#include <IO/Logger.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <IO/LibpllParsers.hpp>

ICCalculator::ICCalculator(const std::string &referenceTreePath,
      const Families &families):
  _rootedReferenceTree(referenceTreePath),
  _referenceTree(_rootedReferenceTree),
  _evaluationTrees()
{
  _readTrees(families);
  _computeQuartets();
  
  _computeScores();  
  printNQuartets(30);
  
}


void ICCalculator::_readTrees(const Families &families)
{
  Logger::timed << "[IC computation] Reading trees..." << std::endl; 
  std::unordered_map<std::string, SPID> speciesLabelToSpid;
  
  for (auto speciesLabel: _rootedReferenceTree.getLabels(true)) {
    auto spid = static_cast<SPID>(speciesLabelToSpid.size());
    speciesLabelToSpid[speciesLabel] = spid;
    _allSPID.insert(spid);
    _spidToString.push_back(speciesLabel);
  }
  _taxaNumber = _allSPID.size();
  _refIdToSPID.resize(_rootedReferenceTree.getLeavesNumber());
  for (auto &leaf: _rootedReferenceTree.getLeaves()) {
    _refIdToSPID[leaf->node_index] = speciesLabelToSpid[std::string(leaf->label)];
  }
  for (auto &family: families) {
    GeneSpeciesMapping mappings;
    mappings.fill(family.mappingFile, family.startingGeneTree);
    auto evaluationTree = std::make_unique<PLLUnrootedTree>(family.startingGeneTree);
    for (auto &leaf: evaluationTree->getLeaves()) {
      leaf->clv_index = speciesLabelToSpid[mappings.getSpecies(std::string(leaf->label))];
    }
    _evaluationTrees.push_back(std::move(evaluationTree));
  }
}


static void getTaxaUnderNode(pll_unode_t *node, TaxaSet &taxa)
{
  if (node->next) {
    getTaxaUnderNode(node->next->back, taxa);
    getTaxaUnderNode(node->next->next->back, taxa);
  } else {
    taxa.insert(node->clv_index);
  }
}

void ICCalculator::_computeQuartets()
{
  Logger::timed << "[IC computation] Computing quartets..." << std::endl; 
  unsigned int quartetsNumber = pow(_allSPID.size(), 4);
  _quartetOccurences = std::vector<SPID>(quartetsNumber, SPID());

  unsigned int printEvery = _evaluationTrees.size() >= 1000 ?
    _evaluationTrees.size() / 10 : 1000;
  if (printEvery < 10) {
    printEvery = 999999999;
  }
  unsigned int plop = 0;
  for (auto &evaluationTree: _evaluationTrees) {
    if (++plop % printEvery == 0) {
      Logger::timed << "    Processed " << plop << "/" << _evaluationTrees.size() << " trees" << std::endl;
    }
    for (auto v: evaluationTree->getInnerNodes()) {
      TaxaSet subtreeTaxa[3];
      getTaxaUnderNode(v->back, subtreeTaxa[0]);
      getTaxaUnderNode(v->next->back, subtreeTaxa[1]);
      getTaxaUnderNode(v->next->next->back, subtreeTaxa[2]);
      for (unsigned int i = 0; i < 3; ++i) {
        for (auto a: subtreeTaxa[i%3]) {
          for (auto b: subtreeTaxa[i%3]) {
            if (a == b) {
              continue;
            }   
            for (auto c: subtreeTaxa[(i+1)%3]) {
              for (auto d: subtreeTaxa[(i+2)%3]) {
                // cover all equivalent quartets
                // we need to swap ab|cd-cd|ab and c-d
                // because of the travserse, we do not need to swap a-b
                _quartetOccurences[_getLookupIndex(a,b,c,d)]++;
                _quartetOccurences[_getLookupIndex(a,b,d,c)]++;
                _quartetOccurences[_getLookupIndex(c,d,a,b)]++;
                _quartetOccurences[_getLookupIndex(d,c,a,b)]++;
              }
            }
          }
        }
      }
    }
  }
}

void ICCalculator::_computeScores()
{
  Logger::timed << "[IC computation] Computing scores..." << std::endl; 
  for (auto u: _referenceTree.getInnerNodes()) {
    for (auto v: _referenceTree.getInnerNodes()) {
      _processNodePair(u, v);
    }
  }
}
  

static SPIDSet _getSPIDClade(
    pll_unode_t *u,
    const std::vector<SPID> &idToSPID)
{
  SPIDSet set;
  for (auto id: PLLUnrootedTree::getClade(u)) {
    set.insert(idToSPID[id]);
  }
  return set;
}

void ICCalculator::_processNodePair(pll_unode_t *u, pll_unode_t *v)
{
  if (u == v) {
    return;
  }
  PLLUnrootedTree::orientTowardEachOther(&u, &v);
  assert(u != v);
  assert(u->next && v->next);
  
  std::array<pll_unode_t*, 4> referenceSubtrees;
  referenceSubtrees[0] = u->next->back;
  referenceSubtrees[1] = u->next->next->back;
  referenceSubtrees[2] = v->next->back;
  referenceSubtrees[3] = v->next->next->back;
  
  std::array<SPIDSet, 4> referenceMetaQuartet;
  for (unsigned int i = 0; i < 4; ++i) {
    referenceMetaQuartet[i] = _getSPIDClade(referenceSubtrees[i], 
        _refIdToSPID);
  }
}

void ICCalculator::printNQuartets(unsigned int n) {
  unsigned int idx = 0;
  for (auto a: _allSPID) {
    for (auto b: _allSPID) {
      if (a == b) {
        continue;
      }
      for (auto c: _allSPID) {
        if (c == a || c == b) {
          continue;
        } 
        for (auto d: _allSPID) {
          if (d == a || d == b || d == c) {
            continue;
          }
          _printQuartet(a, b, c, d);
          
          if (idx++ > n) {
            Logger::info << "Number of species: " << _allSPID.size() << std::endl;
            return;
          }
        }
      }
    }
  }
}

void ICCalculator::_printQuartet(SPID a, SPID b, SPID c, SPID d)
{
  auto astr = _spidToString[a];
  auto bstr = _spidToString[b];
  auto cstr = _spidToString[c];
  auto dstr = _spidToString[d];
  Logger::info << astr << "-" << bstr << " | " << cstr << "-" << dstr;
  //Logger::info << a << "-" << b << " | " << c << "-" << d;
  unsigned int idx[3]; 
  idx[0] = _getLookupIndex(a,b,c,d);
  idx[1] = _getLookupIndex(a,c,b,d);
  idx[2] = _getLookupIndex(a,d,c,b);
  unsigned int occurences[3];
  unsigned int sum = 0;
  for (unsigned int i = 0; i < 3; ++i) {
    occurences[i] = _quartetOccurences[idx[i]];
    sum += occurences[i];
  }
  double frequencies[3];
  for (unsigned int i = 0; i < 3; ++i) {
    frequencies[i] = double(occurences[i]) / double(sum);
    Logger::info << " q" << i << "=" << frequencies[i] << ",";
    //Logger::info << " q" << i << "=" << occurences[i] << ",";
  }
  Logger::info << std::endl;


}




