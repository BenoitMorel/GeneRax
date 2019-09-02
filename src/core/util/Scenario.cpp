#include "util/Scenario.hpp"
#include <IO/Logger.hpp>

const char *Scenario::eventNames[]  = {"S", "SL", "D", "T", "TL", "None", "Invalid"};


void Scenario::addEvent(ReconciliationEventType type, unsigned int geneNode, unsigned int speciesNode) {
  addTransfer(type, geneNode, speciesNode, INVALID);
}
  
void Scenario::addTransfer(ReconciliationEventType type, 
  unsigned int geneNode, 
  unsigned int speciesNode, 
  unsigned int destSpeciesNode)
{
  Event event;
  event.type = type;
  event.geneNode = geneNode;
  event.speciesNode = speciesNode;
  event.destSpeciesNode = destSpeciesNode;
  _events.push_back(event);
  assert(static_cast<int>(type) >= 0);
  _eventsCount[static_cast<unsigned int>(type)] ++;
  if (_geneIdToEvent.size() <= static_cast<size_t>(geneNode)) {
    _geneIdToEvent.resize(geneNode + 1);
  }
  _geneIdToEvent[geneNode] = event;
}

void Scenario::saveEventsCounts(const std::string &filename, bool masterRankOnly) {
  ParallelOfstream os(filename, masterRankOnly);
  for (unsigned int i = 0; i < static_cast<unsigned int>(EVENT_Invalid); ++i) {
    os << eventNames[i] << ":" << _eventsCount[i] << std::endl;
  }
}


void Scenario::recursivelySaveReconciliationsNHX(pll_unode_t *node, ParallelOfstream &os)
{
  if(node->next) {
    os << "(";
    recursivelySaveReconciliationsNHX(node->next->back, os);
    os << ",";
    recursivelySaveReconciliationsNHX(node->next->next->back, os);
    os << ")";
  } 
  if (node->label) {
    os << node->label;
  } else {
    os << "n" << node->node_index; 
  }
  os << ":" << node->length;
  Event event = _geneIdToEvent[node->node_index];
  if (event.speciesNode != INVALID) {
    os << "[&&NHX";
    if (_speciesTree->nodes[event.speciesNode]->label) {
      os << ":S=" << _speciesTree->nodes[event.speciesNode]->label;
    }
    os << ":D=" << (event.type == EVENT_D ? "Y" : "N" );
    os << ":H=" << (event.type == EVENT_T || event.type == EVENT_TL ? "Y" : "N" );
    if (event.type == EVENT_T || event.type == EVENT_TL) {
      assert(_speciesTree->nodes[event.speciesNode]->label);
      assert(_speciesTree->nodes[event.destSpeciesNode]->label);
      os << "@" << _speciesTree->nodes[event.speciesNode]->label;
      os << "@" << _speciesTree->nodes[event.destSpeciesNode]->label;
    }
    os << ":B=" << node->length;
    os << "]";
  }
}
  
void Scenario::saveReconciliationsNHX(const std::string &filename, bool masterRankOnly)
{
  ParallelOfstream os(filename, masterRankOnly);
  os << "(";
  recursivelySaveReconciliationsNHX(_geneRoot, os);
  os << ",";
  recursivelySaveReconciliationsNHX(_geneRoot->back, os);
  os << ");";
}

void Scenario::saveReconciliationsRecPhyloXML(const std::string &filename, bool masterRankOnly)
{
  ParallelOfstream os(filename, masterRankOnly);
  os << "<recPhylo xsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\">" << std::endl;
  
  os << "</recPhylo";
}

void Scenario::saveReconciliations(const std::string &filename, ReconciliationFormat format, bool masterRankOnly)
{
  switch (format) {
  case NHX:
    saveReconciliationsNHX(filename, masterRankOnly);
    break;
  case RecPhyloXML:
    saveReconciliationsRecPhyloXML(filename, masterRankOnly);
    break;
  }
}

