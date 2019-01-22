#include "Scenario.hpp"

void Scenario::addEvent(EventType type, int geneNode, int speciesNode) {
  _os << "(" << eventNames[int(type)] << " " << geneNode << " " << speciesNode << ") ";
  _eventsCount[int(type)] ++;
}

