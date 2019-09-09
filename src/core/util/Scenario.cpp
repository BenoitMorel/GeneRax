#include "util/Scenario.hpp"
#include <IO/Logger.hpp>
#include <IO/ReconciliationWriter.hpp>
#include <IO/ParallelOfstream.hpp>

const char *Scenario::eventNames[]  = {"S", "SL", "D", "T", "TL", "Leaf", "Invalid"};


void Scenario::addEvent(ReconciliationEventType type, unsigned int geneNode, unsigned int speciesNode, unsigned int destSpeciesNode) 
{
  
  addTransfer(type, geneNode, speciesNode, INVALID, destSpeciesNode);
}
  
void Scenario::addTransfer(ReconciliationEventType type, 
  unsigned int geneNode, 
  unsigned int speciesNode, 
    unsigned int transferedGeneNode,
  unsigned int destSpeciesNode)
{
  Event event;
  event.type = type;
  event.geneNode = geneNode;
  event.speciesNode = speciesNode;
  event.transferedGeneNode = transferedGeneNode;
  event.destSpeciesNode = destSpeciesNode;
  _events.push_back(event);
  assert(static_cast<int>(type) >= 0);
  _eventsCount[static_cast<unsigned int>(type)] ++;
  if (_geneIdToEvents.size() <= static_cast<size_t>(geneNode)) {
    _geneIdToEvents.resize(geneNode + 1);
  }
  _geneIdToEvents[geneNode].push_back(event);
}

void Scenario::saveEventsCounts(const std::string &filename, bool masterRankOnly) {
  ParallelOfstream os(filename, masterRankOnly);
  for (unsigned int i = 0; i < static_cast<unsigned int>(EVENT_Invalid); ++i) {
    os << eventNames[i] << ":" << _eventsCount[i] << std::endl;
  }
}

void Scenario::saveReconciliation(const std::string &filename, ReconciliationFormat format, bool masterRankOnly)
{
  switch (format) {
  case NHX:
    ReconciliationWriter::saveReconciliationNHX(_speciesTree, _geneRoot, _virtualRootIndex, _geneIdToEvents, filename, masterRankOnly);
    break;
  case RecPhyloXML:
    ReconciliationWriter::saveReconciliationRecPhyloXML(_speciesTree, _geneRoot, _virtualRootIndex, _geneIdToEvents, filename, masterRankOnly);
    break;
  }
}

