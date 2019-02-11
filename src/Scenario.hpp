#pragma once

#include <vector>
#include <string>
extern "C" {
#include <pll.h>
}
#include <IO/ParallelOfstream.hpp>

using namespace std;

class Scenario {
public:  
  enum EventType {
    S = 0 , SL, D, T, TL, None, Invalid
  };

  struct Event {
    Event(): type(S), geneNode(-1), speciesNode(-1), fromNode(-1), toNode(-1) {}
    EventType type;
    int geneNode;
    int speciesNode;
    int fromNode;
    int toNode; // for transfers
  };


  Scenario(): _eventsCount(int(Invalid), 0), _geneRoot(0) {}
 
  void setGeneRoot(pll_unode_t *geneRoot) {_geneRoot = geneRoot;}
  
  void setSpeciesTree(pll_rtree_t *speciesTree) {_speciesTree = speciesTree;}

  void addEvent(EventType type, int geneNode, int speciesNode);
  void addTransfer(EventType type, int geneNode, int speciesNode, int from, int to);

  void saveEventsCounts(const string &filename); 

  void saveTreeWithEvents(const string &filename);

private:
  static const char *eventNames[];
  vector<Event> _events;
  vector<int> _eventsCount;
  vector<Event> _geneIdToEvent;
  pll_unode_t *_geneRoot;
  pll_rtree_t *_speciesTree;
  void recursivelySaveTreeWithEvents(pll_unode_t *node, ParallelOfstream &os);


};



