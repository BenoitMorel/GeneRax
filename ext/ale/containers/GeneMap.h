// Treerecs – Copyright © by INRIA – All rights reserved – 2017/
// Created by Nicolas Comte on 05/12/16.
//

#ifndef PHYLASOLVER_GENEMAP_H
#define PHYLASOLVER_GENEMAP_H

#include <Bpp/Phyl/Tree/PhyloNode.h>
#include <iostream>
#include <cassert>
#include <memory>
#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <ale/tools/Utils.h>

template<class Species, class Gene>
class GeneMap: protected std::unordered_map<Species, std::list<Gene>>
{
  /*!
   * \class GeneMap
   * \brief Dictionnary associating genes to species.
   * \note We can use std::mulimap<Species, Gene> but we are in a case where multiplicity of items (Genes) per key
   * (Species) is frequent and could be important.
   * Suppose that the number of bins is N and the average number of items per bin is M.
   * A std::multimap<Key, Val> typically uses an RB tree with duplicate keys.
   *    Fetch is O(log N + log M)
   *    Insert is O(log N + log M)
   *    Delete is O(log N + log M)
   *    Iteration is O(1)
   *
   * A std::map<Key, std::vector<Val>> typically uses an RB tree with unique keys.
   *    Fetch is O(log N)
   *    Insert is O(log N)
   *    Delete is O(log N)
   *    Iteration is O(1)
   *
   * The difference between map and multimap is not worth talking about unless M is very large. And that is the point:
   * M can be important.
   */
  std::unordered_map<Gene, Species> geneToSpecies_;

public:
/****************
 * Getters
 */
  /// \brief Generate an unordered_map from the current GeneMap.
  /// \return an std::unordered_map with species as keys and a list of associated genes as values.
  std::unordered_map<Species, std::list<Gene>> u_map() const;

  /// \brief Check and returns (by true or false) if the map gives genes to the given species.
  bool hasGenes(const Species& species) const;

  /// \brief Check if a gene belongs to a given species.
  bool speciesHasGene(const Species& species, const Gene& gene) const;

  /// \brief Returns the number of genes for a given species.
  /// \return Number of genes (std::size_t) associated with a given species.
  ///         returns 0 if their is no gene associated with a species even if the species does not exists in the map.
  std::size_t ngenes(const Species& species) const;

  /// \brief Returns the number of species mapped in the genemap.
  std::size_t nspecies() const;

  /// \brief Returns species of the SpeciesGeneMap into a vector of Species.
  std::vector<Species> getSpecies() const;

  /// \brief Returns a vector of genes associated with a given species.
  std::vector<Gene> getGenes(const Species& species) const;

  /// \brief Returns a vector of mapped genes.
  std::vector<Gene> getGenes(void) const;

  /// \brief Check if a gene has been mapped.
  bool hasGene(const Gene& gene) const;

  /// \brief Returns the species associated with a given gene.
  Species getAssociatedSpecies(const Gene &gene) const;

  bool comp(const GeneMap<Species, Gene> other) const;

/****************
 * Setters
 */
  /// \brief Delete a gene from the map.
  void deleteGene(const Gene& gene);

  /// Delete a gene from the map for a given species.
  void deleteGene(const Species& species, const Gene& gene);

  /// Delete a species (and its associated genes).
  void deleteSpecies(const Species& species);


  /// Add a gene in the map for a given species.
  void addGene(const Species& species, const Gene& gene);

  /// Add a gene in the map for a given species.
  void setGene(const Species& species, const Gene& gene);

  /// Add genes in the map for a given species.
  void addGenes(const Species& species, const std::vector<Gene>& genes);

  /// Add map's informations.
  void addMap(const GeneMap<Species, Gene>& map);
};

template<class Species, class Gene>
std::unordered_map<Species, std::list<Gene>> GeneMap<Species, Gene>::u_map() const {
  return std::unordered_map<Species, std::list<Gene>>(*this);
}

template<class Species, class Gene>
bool GeneMap<Species, Gene>::hasGenes(const Species &species) const {
  if(this->find(species) == this->end())
    return false;

  return(not this->at(species).empty());
}

template<class Species, class Gene>
bool GeneMap<Species, Gene>::speciesHasGene(const Species &species, const Gene &gene) const {
  auto it = geneToSpecies_.find(gene);
  if(it == geneToSpecies_.end())
    return false;
  else
    return geneToSpecies_.at(gene) == species;
}

template<class Species, class Gene>
std::size_t GeneMap<Species, Gene>::ngenes(const Species &species) const {
  if(not this->hasGenes(species))
    return 0;

  return(this->at(species).size());
}

template<class Species, class Gene>
std::size_t GeneMap<Species, Gene>::nspecies() const {
  return this->size();
}

template<class Species, class Gene>
std::vector<Species> GeneMap<Species, Gene>::getSpecies() const {
  std::vector<Species> species;
  species.reserve(this->nspecies());
  for(auto it = this->begin(); it != this->end(); it++){
    species.push_back(it->first);
  }
  return species;
}

template<class Species, class Gene>
std::vector<Gene> GeneMap<Species, Gene>::getGenes(const Species &species) const {
  if(hasGenes(species)) {
    auto& genes = this->at(species);
    return {std::begin(genes), std::end(genes)};
  }
  return {};
}

template<class Species, class Gene>
std::vector<Gene> GeneMap<Species, Gene>::getGenes(void) const {
  std::vector<Gene> genes;
  genes.reserve(geneToSpecies_.size());
  for(auto& pair: geneToSpecies_) genes.push_back(pair.first);
  return genes;
}

template<class Species, class Gene>
bool GeneMap<Species, Gene>::hasGene(const Gene &gene) const {
  for(auto it = this->begin(); it != this->end(); it++){
    if(std::find(it->second.begin(), it->second.end(), gene) != it->second.end())
      return true;
  }
  assert(geneToSpecies_.find(gene) == geneToSpecies_.end());
  return false;
}

template<class Species, class Gene>
Species GeneMap<Species, Gene>::getAssociatedSpecies(const Gene &gene) const {
  if(geneToSpecies_.find(gene) != geneToSpecies_.end()) {
    assert(speciesHasGene(geneToSpecies_.at(gene), gene));
    return geneToSpecies_.at(gene);
  } else {
    std::cerr << "Error, gene " << gene << " does not exists in the map." << std::endl;
    if(hasGene(gene))
      std::cerr << "Warning : the map is not consistent, the gene " << gene
                << " had an association." << std::endl;
    exit(EXIT_FAILURE);
  }
}

template<class Species, class Gene>
bool GeneMap<Species, Gene>::comp(const GeneMap<Species, Gene> other) const {
  if(other.nspecies() != this->nspecies()) return false;

  for(auto this_it = this->begin(); this_it != this->end(); this_it++){
    const auto& species = this_it->first;
    const auto& this_associated_genes_list = this_it->second;
    std::vector<Gene> this_associated_genes{this_associated_genes_list.begin(), this_associated_genes_list.end()};

    if(other.find(species) == other.end()) return false;

    auto other_associated_genes_list = other.getGenes(species);
    std::vector<Gene> other_associated_genes{other_associated_genes_list.begin(), other_associated_genes_list.end()};

    std::sort(this_associated_genes.begin(), this_associated_genes.end());
    std::sort(other_associated_genes.begin(), other_associated_genes.end());

    bool same_associated_genes
        = Utils::comp_all(
            this_associated_genes.begin(), this_associated_genes.end()
            , other_associated_genes.begin(), other_associated_genes.end());

    if(not same_associated_genes) return false;
  }

  return true;
}

template<class Species, class Gene>
void GeneMap<Species, Gene>::deleteGene(const Gene &gene) {
  deleteGene(getAssociatedSpecies(gene), gene);
}

template<class Species, class Gene>
void GeneMap<Species, Gene>::deleteGene(const Species &species, const Gene &gene) {
  if(not this->hasGenes(species)){
    std::cerr << "Error, the species " << species << " is not mapped with the gene " << gene << std::endl;
  }
  auto& genes = this->at(species);
  genes.erase(std::remove(genes.begin(), genes.end(), gene));
  geneToSpecies_.erase(gene);
}

template<class Species, class Gene>
void GeneMap<Species, Gene>::deleteSpecies(const Species &species) {
  for(auto gene: this->getGenes(species)) //First delete all genes
    deleteGene(species, gene);
  this->erase(species); //Delete species then.
}

template<class Species, class Gene>
void GeneMap<Species, Gene>::addGene(const Species &species, const Gene &gene) {
  if(geneToSpecies_.find(gene) != geneToSpecies_.end()){
    //std::cerr << __FUNCTION__ << ": gene " << gene << " already exists for the given species " << getAssociatedSpecies(
    //    gene) << "." << std::endl;
    //assert(false);
    //exit(EXIT_FAILURE);
    deleteGene(getAssociatedSpecies(gene), gene);
  }

  (*this)[species].push_back(gene);
  geneToSpecies_[gene] = species;
}

template<class Species, class Gene>
void GeneMap<Species, Gene>::setGene(const Species &species, const Gene &gene) {
  if(geneToSpecies_.find(gene) != geneToSpecies_.end()){
    deleteGene(getAssociatedSpecies(gene), gene);
  }

  (*this)[species].push_back(gene);
  geneToSpecies_[gene] = species;
}

template<class Species, class Gene>
void GeneMap<Species, Gene>::addGenes(const Species &species, const std::vector<Gene> &genes) {
  for(auto& gene: genes) {
    this->addGene(species, gene);
    geneToSpecies_[gene] = species;
  }
}

template<class Species, class Gene>
void GeneMap<Species, Gene>::addMap(const GeneMap<Species, Gene> &map) {
  for(auto& pair: map){
    for(auto& gene: pair.second){
      if(not this->speciesHasGene(pair.first, gene))
        this->addGene(pair.first, gene);
    }
  }
}

template <typename Species, typename Gene>
inline bool operator==(const GeneMap<Species, Gene>& map0, const GeneMap<Species, Gene>& map1)
{ return map0.comp(map1); }

using SpeciesGeneMap = GeneMap<std::shared_ptr<bpp::PhyloNode>, std::shared_ptr<bpp::PhyloNode>> ;

#endif //PHYLASOLVER_GENEMAP_H
