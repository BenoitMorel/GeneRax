#include "Clade.hpp"
#include <algorithm>
#include <functional>
#include <trees/PLLUnrootedTree.hpp>
#include <trees/PLLRootedTree.hpp>
#include <IO/GeneSpeciesMapping.hpp>

static unsigned int hashTwoValues(unsigned int a, unsigned int b)
{
  /*
  std::hash<size_t> hash_fn;
  return hash_fn(a*106039+b);
  */
  return a*106039+b;
}

Clade::Clade():
  _hash(0)
{

}

void Clade::mergeWith(const Clade &clade)
{
  auto size = _ids.size();
  _ids.insert(clade._ids.begin(), clade._ids.end());
  if (size != _ids.size()) {
    _recomputeHash();
  }
}

void Clade::addId(unsigned int id)
{
  auto size = _ids.size();
  _ids.insert(id);
  if (size != _ids.size()) {
    _recomputeHash();
  }
}

void Clade::_recomputeHash()
{
  _hash = 42;
  for (auto id: _ids) {
    _hash = hashTwoValues(_hash, id);
  }
}

Clade Clade::getComplement(Clade &allTaxa)
{
  Clade complement;
  std::set_difference(allTaxa._ids.begin(), allTaxa._ids.end(),
      _ids.begin(), _ids.end(), 
      std::inserter(complement._ids, complement._ids.end()));
  complement._recomputeHash();
  return complement;
}
  
Clade Clade::getMaximumClade(PLLRootedTree &tree)
{
  Clade clade;
  for (auto &p: tree.getLabelToIntMap()) {
    clade._ids.insert(p.second);
  }
  clade._recomputeHash();
  return clade;
}


CladeSet Clade::buildCladeSet(PLLRootedTree &tree)
{
  auto postOrderNodes = tree.getPostOrderNodes();
  std::vector<Clade> clades(postOrderNodes.size());
  for (auto node: postOrderNodes) {
    auto &clade = clades[node->node_index];
    if (!node->left) {
      // leaf node
      clade.addId(node->node_index);
    } else {
      // internal node
      const auto &cladeLeft = clades[node->left->node_index];
      const auto &cladeRight = clades[node->right->node_index];
      clade.mergeWith(cladeLeft);
      clade.mergeWith(cladeRight);
    }
  }
  auto maximumClade = clades.back();
  auto emptyClade = Clade();
  CladeSet res;
  for (auto &clade: clades) {
    res.insert(clade.getHash());
  }
  res.erase(maximumClade.getHash());
  res.erase(emptyClade.getHash());
  return res;
}


CladeSet Clade::buildCladeSet(PLLUnrootedTree &tree,
      const GeneSpeciesMapping &geneSpeciesMapping,
      const StringToUintMap &speciesLabelToInt)
{
  auto postOrderNodes = tree.getPostOrderNodes();
  std::vector<Clade> clades(postOrderNodes.size());
  for (auto node: postOrderNodes) {
    auto &clade = clades[node->node_index];
    if (!node->next) {
      // leaf node
      auto id = speciesLabelToInt.at(geneSpeciesMapping.getSpecies(node->label));
      clade.addId(id);
    } else {
      // internal node
      const auto &cladeLeft = clades[node->next->back->node_index];
      const auto &cladeRight = clades[node->next->next->back->node_index];
      clade.mergeWith(cladeLeft);
      clade.mergeWith(cladeRight);
    }
  }
  bool bipartitionsOnly = true;
  CladeSet res;
  if (bipartitionsOnly) {
    for (auto node: postOrderNodes) {
      auto &c1 = clades[node->node_index];
      auto &c2 = clades[node->back->node_index];
      if (intersection_size(c1._ids, c2._ids) == 0) {
        res.insert(c1.getHash());  
        res.insert(c2.getHash());  
      }
    }
  } else {
    for (auto &clade: clades) {
      res.insert(clade.getHash());
    }
  }
  return res;
}
