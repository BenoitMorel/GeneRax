#pragma once

#include <util/Scenario.hpp>
#include <util/enums.hpp>
#include <vector>
#include <string>
typedef struct pll_unode_s pll_unode_t;
typedef struct pll_rtree_s pll_rtree_t;


class ReconciliationWriter {
public:
  static void saveReconciliationNHX(pll_rtree_t *speciesTree,  
      pll_unode_t *geneRoot, 
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      const std::string &filename,
      bool masterRankOnly); 

  static void saveReconciliationRecPhyloXML(pll_rtree_t *speciesTree,  
      pll_unode_t *geneRoot, 
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      const std::string &filename,
      bool masterRankOnly); 
};




