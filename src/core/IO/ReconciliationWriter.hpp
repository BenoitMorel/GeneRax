#pragma once

#include <util/Scenario.hpp>
#include <util/enums.hpp>
#include <vector>
#include <string>
typedef struct corax_unode_s corax_unode_t;
typedef struct corax_rtree_s corax_rtree_t;
class ParallelOfstream;

class ReconciliationWriter {
public:
  ReconciliationWriter() = delete;

  static void saveReconciliationNHX(corax_rtree_t *speciesTree,  
      corax_unode_t *geneRoot, 
      unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      ParallelOfstream &os);

  static void saveReconciliationRecPhyloXML(corax_rtree_t *speciesTree,  
      unsigned int geneNode, 
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      ParallelOfstream &os);
  
  static void saveReconciliationNewickEvents(corax_unode_t *geneRoot, 
      unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      ParallelOfstream &os);
};




