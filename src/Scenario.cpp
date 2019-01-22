#include "Scenario.hpp"

const char *Scenario::eventNames[]  = {"S", "SL", "D", "T", "None", "Invalid"};


void Scenario::addEvent(EventType type, int geneNode, int speciesNode) {
  Event event;
  event.type = type;
  event.geneNode = geneNode;
  event.speciesNode = speciesNode;
  _events.push_back(event);
  _eventsCount[int(type)] ++;
}

