#include "Scenario.hpp"
#include <IO/Logger.hpp>

const char *Scenario::eventNames[]  = {"S", "SL", "D", "T", "TL", "None", "Invalid"};


void Scenario::addEvent(EventType type, int geneNode, int speciesNode) {
  addTransfer(type, geneNode, speciesNode, -1, -1);
}
  
void Scenario::addTransfer(EventType type, int geneNode, int speciesNode, int from, int to)
{
  Event event;
  event.type = type;
  event.geneNode = geneNode;
  event.speciesNode = speciesNode;
  event.fromNode = from;
  event.toNode = to;
  _events.push_back(event);
  _eventsCount[int(type)] ++;
  if (_geneIdToEvent.size() <= (size_t)geneNode) {
    _geneIdToEvent.resize(geneNode + 1);
  }
  _geneIdToEvent[geneNode] = event;
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
  } 
  if (node->label) {
    os << node->label;
  } else {
    os << "n" << node->node_index; 
  }
  Event event = _geneIdToEvent[node->node_index];
  if (event.speciesNode >= 0) {
    os << "[&&NHX";
    if (_speciesTree->nodes[event.speciesNode]->label) {
      os << ":S=" << _speciesTree->nodes[event.speciesNode]->label;
    }
    os << ":D=" << (event.type == D ? "Y" : "N" );
    os << ":H=" << (event.type == T || event.type == TL ? "Y" : "N" );
    if (event.type == T || event.type == TL) {
      os << "@" << _speciesTree->nodes[event.fromNode]->label;
      os << "@" << _speciesTree->nodes[event.toNode]->label;
    }
    os << ":B=" << node->length;
    os << "]";
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

