#pragma once

#include <vector>
#include <string>
#include <util/enums.hpp>
extern "C" {
#include <pll.h>
}

typedef struct pll_utree_s pll_utree_t;
typedef struct pll_unode_s pll_unode_t;
typedef struct pll_rtree_s pll_rtree_t;
typedef struct pll_rnode_s pll_rnode_t;

/*
 * Store the set of events that reconciles a gene tree with a species tree
 */
class Scenario {
public:
  static const unsigned int INVALID = static_cast<unsigned int>(-1);

  struct Event {
    Event(): type(EVENT_S), geneNode(INVALID), speciesNode(INVALID), transferedGeneNode(INVALID), destSpeciesNode(INVALID) {}
    ReconciliationEventType type;
    unsigned int geneNode;
    unsigned int speciesNode;
    unsigned int transferedGeneNode; 
    unsigned int destSpeciesNode; // for transfers
    bool isValid() const { return speciesNode != INVALID; }
  };


  Scenario(): _eventsCount(static_cast<unsigned int>(EVENT_Invalid), 0), _geneRoot(0), _virtualRootIndex(-1) {}
 
  void setGeneRoot(pll_unode_t *geneRoot) {_geneRoot = geneRoot;}
  
  void setSpeciesTree(pll_rtree_t *speciesTree) {_speciesTree = speciesTree;}
  void setVirtualRootIndex(unsigned int virtualRootIndex) {_virtualRootIndex = virtualRootIndex;}

  void addEvent(ReconciliationEventType type, unsigned int geneNode, unsigned int speciesNode, unsigned int destSpeciesNode = INVALID);
  void addTransfer(ReconciliationEventType type, 
    unsigned int geneNode, 
    unsigned int speciesNode, 
    unsigned int transferedGeneNode,
    unsigned int destSpeciesNode);

  void saveEventsCounts(const std::string &filename, bool masterRankOnly = true); 

  void saveReconciliation(const std::string &filename, ReconciliationFormat format, bool masterRankOnly = true);

private:
  static const char *eventNames[];
  std::vector<Event> _events;
  std::vector<unsigned int> _eventsCount;
  std::vector<std::vector<Event> > _geneIdToEvents;
  pll_unode_t *_geneRoot;
  pll_rtree_t *_speciesTree;
  unsigned int _virtualRootIndex;
};



