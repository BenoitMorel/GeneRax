#include "ReconciliationWriter.hpp"

#include <IO/ParallelOfstream.hpp>
#include <corax/corax_common.h>
#include <trees/PLLRootedTree.hpp>

static void printEvent(const Scenario::Event &event, pll_rtree_t *speciesTree, pll_unode_t *node, ParallelOfstream &os)
{
  if (event.isValid()) {
    os << "[&&NHX";
    if (speciesTree->nodes[event.speciesNode]->label) {
      os << ":S=" << speciesTree->nodes[event.speciesNode]->label;
    }
    os << ":D=" << (event.type == ReconciliationEventType::EVENT_D ? "Y" : "N" );
    os << ":H=" << 
      (event.type == ReconciliationEventType::EVENT_T || event.type == ReconciliationEventType::EVENT_TL ?
       "Y" : "N" );
    if (event.type == ReconciliationEventType::EVENT_T || event.type == ReconciliationEventType::EVENT_TL) {
      assert(speciesTree->nodes[event.speciesNode]->label);
      assert(speciesTree->nodes[event.destSpeciesNode]->label);
      os << "@" << speciesTree->nodes[event.speciesNode]->label;
      os << "@" << speciesTree->nodes[event.destSpeciesNode]->label;
    }
    os << ":B=" << node->length;
    os << "]";
  }
}

static void recursivelySaveReconciliationsNHX(pll_rtree_t *speciesTree, 
    pll_unode_t *node, 
    bool isVirtualRoot, 
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  
  if(node->next) {
    pll_unode_t *left = nullptr;
    pll_unode_t *right = nullptr;
    if (isVirtualRoot) {
      left = node->next;
      right = node->next->back;
    } else {
      left = node->next->back;
      right = node->next->next->back;
    }
    os << "(";
    recursivelySaveReconciliationsNHX(speciesTree, left, false, geneToEvents, os);
    os << ",";
    recursivelySaveReconciliationsNHX(speciesTree, right, false, geneToEvents, os);
    os << ")";
  } 
  if (node->label) {
    os << node->label;
  } else {
    os << "n" << node->node_index; 
  }
  if (!isVirtualRoot) {
    os << ":" << node->length;
  }
  printEvent(geneToEvents[node->node_index].back(), speciesTree, node, os);
}
  
void ReconciliationWriter::saveReconciliationNHX(pll_rtree_t *speciesTree, 
    pll_unode_t *geneRoot, 
    unsigned int virtualRootIndex,
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os) 
{
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  recursivelySaveReconciliationsNHX(speciesTree, &virtualRoot, true, geneToEvents, os);
  os << ";";
}


static void recursivelySaveSpeciesTreeRecPhyloXML(pll_rnode_t *node, std::string &indent, ParallelOfstream &os)
{
  if (!node) {
    return;
  }
  os << indent << "<clade>" << std::endl;
  indent += "\t";
  os << indent << "\t<name>" << node->label << "</name>" << std::endl;
  recursivelySaveSpeciesTreeRecPhyloXML(node->left, indent, os);
  recursivelySaveSpeciesTreeRecPhyloXML(node->right, indent, os);
  indent.pop_back();
  os << indent << "</clade>" << std::endl;
}

static void saveSpeciesTreeRecPhyloXML(pll_rtree_t *speciesTree, ParallelOfstream &os)
{
  os << "<spTree>" << std::endl;
  os << "<phylogeny>" << std::endl;
  std::string indent = "";
  
  recursivelySaveSpeciesTreeRecPhyloXML(speciesTree->root, indent, os);
  os << "</phylogeny>" << std::endl;
  os << "</spTree>" << std::endl;
}

static void writeEventRecPhyloXML(pll_unode_t *geneTree,
    pll_rtree_t *speciesTree, 
    Scenario::Event  &event,
    const Scenario::Event *previousEvent,
    std::string &indent, 
    ParallelOfstream &os)
{
  auto species = speciesTree->nodes[event.speciesNode];
  pll_rnode_t *speciesOut = 0;
  os << indent << "<eventsRec>" << std::endl;
  bool previousWasTransfer = previousEvent->type 
    == ReconciliationEventType::EVENT_T || previousEvent->type == ReconciliationEventType::EVENT_TL;
  if (previousWasTransfer && geneTree->node_index == previousEvent->rightGeneIndex && event.type 
      != ReconciliationEventType::EVENT_L) {
    auto previousEventSpeciesOut = speciesTree->nodes[previousEvent->destSpeciesNode];
    os << indent << "\t<transferBack destinationSpecies=\"" << previousEventSpeciesOut->label << "\"/>" << std::endl;
  }
  
  switch(event.type) {
  case ReconciliationEventType::EVENT_None:
    assert(geneTree->next == 0);
    assert(species->left == 0 && species->right == 0);
    os << indent << "\t<leaf speciesLocation=\"" << species->label << "\"/>" <<  std::endl;
    break;
  case ReconciliationEventType::EVENT_S: case ReconciliationEventType::EVENT_SL:
    os << indent << "\t<speciation speciesLocation=\"" << species->label << "\"/>" << std::endl;
    break;
  case ReconciliationEventType::EVENT_D:
    os << indent << "\t<duplication speciesLocation=\"" << species->label << "\"/>" << std::endl;
    break;
  case ReconciliationEventType::EVENT_T: case ReconciliationEventType::EVENT_TL:
    speciesOut = speciesTree->nodes[event.speciesNode];
    os << indent << "\t<branchingOut speciesLocation=\"" << speciesOut->label << "\"/>" << std::endl;
    break; 
  case ReconciliationEventType::EVENT_L:
    speciesOut = speciesTree->nodes[event.speciesNode];
    os << indent << "\t<loss speciesLocation=\"" << speciesOut->label << "\"/>" << std::endl;
    break;
  default:
    break;
  }
  os << indent << "</eventsRec>" << std::endl;
}

static void recursivelySaveGeneTreeRecPhyloXML(pll_unode_t *geneTree, 
    bool isVirtualRoot,
    pll_rtree_t *speciesTree, 
    std::vector<std::vector<Scenario::Event> > &geneToEvents,
    const Scenario::Event *previousEvent,
    std::string &indent,
    ParallelOfstream &os)
{
  if (!geneTree) {
    return;
  }
  auto &events = geneToEvents[geneTree->node_index];
  for (unsigned int i = 0; i < events.size() - 1; ++i) {
    os << indent << "<clade>" << std::endl;
    indent += "\t";
    auto &event = events[i];
    os << indent << "<name>" << (geneTree->label ? geneTree->label : "NULL") << "</name>" << std::endl;
    writeEventRecPhyloXML(geneTree, speciesTree, event, previousEvent, indent, os);  
    previousEvent = &event;
    if (event.type == ReconciliationEventType::EVENT_SL || event.type == ReconciliationEventType::EVENT_TL) {
      Scenario::Event loss;
      loss.type = ReconciliationEventType::EVENT_L;
      if (event.type == ReconciliationEventType::EVENT_SL) {
        auto parentSpecies = speciesTree->nodes[event.speciesNode];
        auto lostSpecies = (parentSpecies->left->node_index == event.destSpeciesNode) 
          ? parentSpecies->right : parentSpecies->left;
        loss.speciesNode = lostSpecies->node_index;
      } else if (event.type == ReconciliationEventType::EVENT_TL) {
        loss.speciesNode = event.speciesNode;
      }
      indent += "\t";
      os << indent << "<clade>" << std::endl;
      os << indent << "<name>loss</name>" << std::endl;
      writeEventRecPhyloXML(geneTree, speciesTree, loss, previousEvent, indent, os);
      indent.pop_back();
      os << indent << "</clade>" << std::endl;
    } else {
      assert(false); 
    }
  }

  os << indent << "<clade>" << std::endl;
  indent += "\t";
  Scenario::Event &event = geneToEvents[geneTree->node_index].back();
  os << indent << "<name>" << (geneTree->label ? geneTree->label : "NULL") << "</name>" << std::endl;
  writeEventRecPhyloXML(geneTree, speciesTree, event, previousEvent, indent, os);  


  if (geneTree->next) {
    pll_unode_t *left = nullptr;
    pll_unode_t *right = nullptr;
    if (isVirtualRoot) {
      
      assert(geneTree->next);
      assert(geneTree->next->back);
      left = geneTree->next;
      right = geneTree->next->back;
    } else {
      left = geneTree->next->back;
      right = geneTree->next->next->back;
    }
    
    recursivelySaveGeneTreeRecPhyloXML(left, false, speciesTree, geneToEvents, &event, indent, os);
    recursivelySaveGeneTreeRecPhyloXML(right, false, speciesTree, geneToEvents, &event, indent, os);
  }
  for (unsigned int i = 0; i < events.size() - 1; ++i) {
    indent.pop_back();
    os << indent << "</clade>" << std::endl;
  }
  indent.pop_back();
  os << indent << "</clade>" << std::endl;
}

static void saveGeneTreeRecPhyloXML(pll_unode_t *geneTree,
    unsigned int virtualRootIndex,
    pll_rtree_t *speciesTree,
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  os << "<recGeneTree>" << std::endl;
  os << "<phylogeny rooted=\"true\">" << std::endl;
  std::string indent;
  Scenario::Event noEvent;
  noEvent.type = ReconciliationEventType::EVENT_None;
  pll_unode_t virtualRoot;
  virtualRoot.next = geneTree;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = 0;
  recursivelySaveGeneTreeRecPhyloXML(&virtualRoot, true, speciesTree, geneToEvents, &noEvent, indent, os); 
  os << "</phylogeny>" << std::endl;
  os << "</recGeneTree>" << std::endl;
}

void ReconciliationWriter::saveReconciliationRecPhyloXML(pll_rtree_t *speciesTree, 
    pll_unode_t *geneRoot, 
    unsigned int virtualRootIndex,
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  os << "<recPhylo " << std::endl;
  os << "\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << std::endl;
  os << "\txsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\"" << std::endl;
  os << "\txmlns=\"http://www.recg.org\">" << std::endl;
  saveSpeciesTreeRecPhyloXML(speciesTree, os);
  saveGeneTreeRecPhyloXML(geneRoot, virtualRootIndex, speciesTree, geneToEvents, os);
  os << "</recPhylo>";

}

static void recursivelySaveReconciliationsNewickEvents(
    pll_unode_t *node, 
    bool isVirtualRoot, 
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  
  if(node->next) {
    auto left = node->next->back;
    auto right = node->next->next->back;
    if (isVirtualRoot) {
      left = node->next;
      right = node->next->back;
    }
    os << "(";
    recursivelySaveReconciliationsNewickEvents(
        left, 
        false, 
        geneToEvents, 
        os);
    os << ",";
    recursivelySaveReconciliationsNewickEvents(
        right, 
        false, 
        geneToEvents, 
        os);
    os << ")";
  } 
  if (!node->next) {
    os << node->label;
  } else {
    os << Enums::getEventName(geneToEvents[node->node_index].back().type); 
  }
  if (!isVirtualRoot) {
    os << ":" << node->length;
  }
}
  
void ReconciliationWriter::saveReconciliationNewickEvents(pll_unode_t *geneRoot, 
      unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      ParallelOfstream &os)
{
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  recursivelySaveReconciliationsNewickEvents(&virtualRoot, 
      true, 
      geneToEvent, 
      os);
  os << ";";

}
