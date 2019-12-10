#pragma once

#include <util/Scenario.hpp>
#include <util/enums.hpp>
#include <vector>
#include <string>
typedef struct pll_unode_s pll_unode_t;
typedef struct pll_rtree_s pll_rtree_t;
class ParallelOfstream;

class ReconciliationWriter {
public:
  ReconciliationWriter() = delete;

  static void saveReconciliationNHX(pll_rtree_t *speciesTree,  
      pll_unode_t *geneRoot, 
      unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      ParallelOfstream &os);

  static void saveReconciliationRecPhyloXML(pll_rtree_t *speciesTree,  
      pll_unode_t *geneRoot, 
      unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      ParallelOfstream &os);
};




