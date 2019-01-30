#include "Scenario.hpp"
#include <IO/Logger.hpp>

const char *Scenario::eventNames[]  = {"S", "SL", "D", "T", "TL", "None", "Invalid"};


void Scenario::addEvent(EventType type, int geneNode, int speciesNode) {
  Event event;
  event.type = type;
  event.geneNode = geneNode;
  event.speciesNode = speciesNode;
  _events.push_back(event);
  _eventsCount[int(type)] ++;
  if (_geneIdToEvent.size() <= (size_t)geneNode) {
    _geneIdToEvent.resize(geneNode + 1);
  }
  _geneIdToEvent[geneNode] = type;
}
  
void Scenario::saveEventsCounts(const string &filename) {
  ParallelOfstream os(filename);
  for (int i = 0; i < int(Invalid); ++i) {
    os << eventNames[i] << ":" << _eventsCount[i] << endl;
  }
}


void Scenario::recursivelySaveTreeWithEvents(pll_unode_t *node, ParallelOfstream &os)
{
  if(node->next) {
    os << "(";
    recursivelySaveTreeWithEvents(node->next->back, os);
    os << ",";
    recursivelySaveTreeWithEvents(node->next->next->back, os);
    os << ")";
  } else if (node->label) {
    os << node->label;
  }
  EventType event = _geneIdToEvent[node->node_index];
  if (event == T || event == D || event == TL) {
    os << eventNames[int(event)];
  }

}

void Scenario::saveTreeWithEvents(const string &filename)
{
  ParallelOfstream os(filename);
  os << "(";
  recursivelySaveTreeWithEvents(_geneRoot, os);
  os << ",";
  recursivelySaveTreeWithEvents(_geneRoot->back, os);
  os << ");";
}

