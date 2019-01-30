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
    EventType type;
    int geneNode;
    int speciesNode;
  };


  Scenario(): _eventsCount(int(Invalid), 0), _geneRoot(0) {}
 
  void setGeneRoot(pll_unode_t *geneRoot) {_geneRoot = geneRoot;}

  void addEvent(EventType type, int geneNode, int speciesNode);

  void saveEventsCounts(const string &filename); 

  void saveTreeWithEvents(const string &filename);

private:
  static const char *eventNames[];
  vector<Event> _events;
  vector<int> _eventsCount;
  vector<EventType> _geneIdToEvent;
  pll_unode_t *_geneRoot;

  void recursivelySaveTreeWithEvents(pll_unode_t *node, ParallelOfstream &os);


};



