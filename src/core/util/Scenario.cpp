#include "util/Scenario.hpp"
#include <IO/Logger.hpp>
#include <IO/ReconciliationWriter.hpp>
#include <IO/ParallelOfstream.hpp>
#include <vector>
#include <string>

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

void Scenario::savePerSpeciesEventsCounts(const std::string &filename, bool masterRankOnly) {
  
  ParallelOfstream os(filename, masterRankOnly);
  std::map<std::string, std::vector<unsigned int> > speciesToEventCount;
  std::vector<unsigned int> defaultCount(static_cast<unsigned int>(4), 0);
  for (unsigned int e = 0; e < _speciesTree->tip_count + _speciesTree->inner_count; ++e) {
    assert(_speciesTree->nodes[e]->label);
    speciesToEventCount.insert(std::pair<std::string, std::vector<unsigned int> > (std::string(_speciesTree->nodes[e]->label), defaultCount));
  }
  for (auto &event: _events) {
    auto &eventCount = speciesToEventCount[_speciesTree->nodes[event.speciesNode]->label];
    switch (event.type) {
      case EVENT_S: 
      case EVENT_None:
        eventCount[0]++;
        break;
      case EVENT_SL: 
        eventCount[0]++;
        eventCount[2]++;
      case EVENT_D:
        eventCount[1]++;
        break;
      case EVENT_T:
        eventCount[3]++;
        break;
      case EVENT_TL:
        eventCount[2]++;
        eventCount[3]++;
        break;
      default:
        break;
    }
  }
  for (auto &it: speciesToEventCount) {
    os << it.first << " ";
    for (auto v: it.second) {
      os << v << " ";
    }
    os << "\n";
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
  
void Scenario::saveTransfers(const std::string &filename, bool masterRankOnly)
{
  ParallelOfstream os(filename, masterRankOnly);
  for (auto &event: _events) {
    if (event.type == EVENT_T || event.type == EVENT_TL) {
      os << _speciesTree->nodes[event.speciesNode]->label << " " << _speciesTree->nodes[event.destSpeciesNode]->label << std::endl;
    }
  }
}





