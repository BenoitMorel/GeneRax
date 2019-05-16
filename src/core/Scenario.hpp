#pragma once

#include <vector>
#include <string>
#include <IO/ParallelOfstream.hpp>
extern "C" {
#include <pll.h>
}

typedef struct pll_utree_s pll_utree_t;
typedef struct pll_unode_s pll_unode_t;
typedef struct pll_rtree_s pll_rtree_t;
typedef struct pll_rnode_s pll_rnode_t;

class Scenario {
public:
  static const unsigned int INVALID = static_cast<unsigned int>(-1);

  enum EventType {
    S = 0 , SL, D, T, TL, None, Invalid
  };

  struct Event {
    Event(): type(S), geneNode(INVALID), speciesNode(INVALID), fromNode(INVALID), toNode(INVALID) {}
    EventType type;
    unsigned int geneNode;
    unsigned int speciesNode;
    unsigned int fromNode;
    unsigned int toNode; // for transfers
  };


  Scenario(): _eventsCount(static_cast<unsigned int>(Invalid), 0), _geneRoot(0) {}
 
  void setGeneRoot(pll_unode_t *geneRoot) {_geneRoot = geneRoot;}
  
  void setSpeciesTree(pll_rtree_t *speciesTree) {_speciesTree = speciesTree;}

  void addEvent(EventType type, unsigned int geneNode, unsigned int speciesNode);
  void addTransfer(EventType type, 
    unsigned int geneNode, 
    unsigned int speciesNode, 
    unsigned int from, 
    unsigned int to);

  void saveEventsCounts(const std::string &filename, bool masterRankOnly = true); 

  void saveTreeWithEvents(const std::string &filename, bool masterRankOnly = true);

private:
  static const char *eventNames[];
  std::vector<Event> _events;
  std::vector<unsigned int> _eventsCount;
  std::vector<Event> _geneIdToEvent;
  pll_unode_t *_geneRoot;
  pll_rtree_t *_speciesTree;
  void recursivelySaveTreeWithEvents(pll_unode_t *node, ParallelOfstream &os);


};



