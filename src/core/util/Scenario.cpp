#include "util/Scenario.hpp"
#include <IO/Logger.hpp>
#include <IO/ReconciliationWriter.hpp>
#include <IO/ParallelOfstream.hpp>
#include <vector>
#include <string>
#include <map>


const char *Scenario::eventNames[]  = {"S", "SL", "D", "T", "TL", "L", "Leaf", "Invalid"};


void Scenario::addEvent(ReconciliationEventType type, 
    unsigned int geneNode, 
    unsigned int speciesNode, 
    unsigned int destSpeciesNode) 
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
  addEvent(event);
}
  
void Scenario::addEvent(const Event &event)
{
  _events.push_back(event);
  assert(static_cast<int>(event.type) >= 0);
  _eventsCount[static_cast<unsigned int>(event.type)] ++;
  if (_geneIdToEvents.size() <= static_cast<size_t>(event.geneNode)) {
    _geneIdToEvents.resize(event.geneNode + 1);
  }
  _geneIdToEvents[event.geneNode].push_back(event);

}

void Scenario::saveEventsCounts(const std::string &filename, bool masterRankOnly) {
  ParallelOfstream os(filename, masterRankOnly);
  for (unsigned int i = 0; i < static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid); ++i) {
    os << eventNames[i] << ":" << _eventsCount[i] << std::endl;
  }
}

void Scenario::savePerSpeciesEventsCounts(const std::string &filename, bool masterRankOnly) {
  
  ParallelOfstream os(filename, masterRankOnly);
  std::map<std::string, std::vector<unsigned int> > speciesToEventCount;
  std::vector<unsigned int> defaultCount(static_cast<unsigned int>(4), 0);
  for (unsigned int e = 0; e < _speciesTree->tip_count + _speciesTree->inner_count; ++e) {
    assert(_speciesTree->nodes[e]->label);
    speciesToEventCount.insert({std::string(_speciesTree->nodes[e]->label), defaultCount});
  }
  for (auto &event: _events) {
    auto &eventCount = speciesToEventCount[_speciesTree->nodes[event.speciesNode]->label];
    switch (event.type) {
      case ReconciliationEventType::EVENT_S: 
      case ReconciliationEventType::EVENT_None:
        eventCount[0]++;
        break;
      case ReconciliationEventType::EVENT_SL: 
        eventCount[0]++;
        eventCount[2]++;
        break;
      case ReconciliationEventType::EVENT_D:
        eventCount[1]++;
        break;
      case ReconciliationEventType::EVENT_T:
        eventCount[3]++;
        break;
      case ReconciliationEventType::EVENT_TL:
        eventCount[2]++;
        eventCount[3]++;
        break;
      case ReconciliationEventType::EVENT_L: case ReconciliationEventType::EVENT_Invalid:
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
  ParallelOfstream os(filename, masterRankOnly);
  saveReconciliation(os, format);
}

void Scenario::saveReconciliation(ParallelOfstream &os, ReconciliationFormat format)
{
  switch (format) {
  case ReconciliationFormat::NHX:
    ReconciliationWriter::saveReconciliationNHX(_speciesTree, 
        _geneRoot, 
        _virtualRootIndex, 
        _geneIdToEvents, 
        os);
    break;
  case ReconciliationFormat::RecPhyloXML:
    ReconciliationWriter::saveReconciliationRecPhyloXML(_speciesTree, 
        _geneRoot, 
        _virtualRootIndex, 
        _geneIdToEvents, 
        os);
    break;
  }
}
  
void Scenario::saveTransfers(const std::string &filename, bool masterRankOnly)
{
  ParallelOfstream os(filename, masterRankOnly);
  for (auto &event: _events) {
    if (event.type == ReconciliationEventType::EVENT_T || event.type == ReconciliationEventType::EVENT_TL) {
      os << _speciesTree->nodes[event.speciesNode]->label << " " 
        << _speciesTree->nodes[event.destSpeciesNode]->label << std::endl;
    }
  }
}
void Scenario::initBlackList(unsigned int genesNumber, unsigned int speciesNumber)
{
  std::vector<int>  blackListAux(speciesNumber, false);
  _blacklist = std::make_unique<ScenarioBlackList>(genesNumber, blackListAux);
  
}

void Scenario::blackList(unsigned int geneNode, unsigned int speciesNode)
{
  if (_blacklist && _blacklist->size() > geneNode) { // not true for virtual nodes
    (*_blacklist)[geneNode][speciesNode] = true;
  }
}
bool Scenario::isBlacklisted(unsigned int geneNode, unsigned int speciesNode)
{
  if (_blacklist && _blacklist->size() > geneNode) { // not true for virtual nodes
    return (*_blacklist)[geneNode][speciesNode];
  }
  return false;
}

void Scenario::resetBlackList()
{
  for (auto &aux: *_blacklist) {
    for (unsigned int i = 0; i < aux.size(); ++i) {
      aux[i] = false;
    }
  }
}
