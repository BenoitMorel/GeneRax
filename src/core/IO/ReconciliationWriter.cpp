#include "ReconciliationWriter.hpp"

#include <IO/ParallelOfstream.hpp>
#include <corax/corax.h>
#include <trees/PLLRootedTree.hpp>

static void printEvent(const Scenario::Event &event, corax_rtree_t *speciesTree, corax_unode_t *node, ParallelOfstream &os)
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

static void recursivelySaveReconciliationsNHX(corax_rtree_t *speciesTree, 
    corax_unode_t *node, 
    bool isVirtualRoot, 
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  
  if(node->next) {
    corax_unode_t *left = nullptr;
    corax_unode_t *right = nullptr;
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
  
void ReconciliationWriter::saveReconciliationNHX(corax_rtree_t *speciesTree, 
    corax_unode_t *geneRoot, 
    unsigned int virtualRootIndex,
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os) 
{
  corax_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  recursivelySaveReconciliationsNHX(speciesTree, &virtualRoot, true, geneToEvents, os);
  os << ";";
}


static void recursivelySaveSpeciesTreeRecPhyloXML(corax_rnode_t *node, std::string &indent, ParallelOfstream &os)
{
  if (!node) {
    return;
  }
  os << indent << "<clade>" << std::endl;
  indent += "\t";
  std::string label(node->label);
  if (!label.size()) {
    label = "NULL";
  }
  os << indent << "\t<name>" << label << "</name>" << std::endl;
  recursivelySaveSpeciesTreeRecPhyloXML(node->left, indent, os);
  recursivelySaveSpeciesTreeRecPhyloXML(node->right, indent, os);
  indent.pop_back();
  os << indent << "</clade>" << std::endl;
}

static void saveSpeciesTreeRecPhyloXML(corax_rtree_t *speciesTree, ParallelOfstream &os)
{
  os << "<spTree>" << std::endl;
  os << "<phylogeny>" << std::endl;
  std::string indent = "";
  
  recursivelySaveSpeciesTreeRecPhyloXML(speciesTree->root, indent, os);
  os << "</phylogeny>" << std::endl;
  os << "</spTree>" << std::endl;
}

static void writeEventRecPhyloXML(unsigned int geneIndex,
    corax_rtree_t *speciesTree, 
    Scenario::Event  &event,
    const Scenario::Event *previousEvent,
    std::string &indent, 
    ParallelOfstream &os)
{
  auto species = speciesTree->nodes[event.speciesNode];
  corax_rnode_t *speciesOut = 0;
  os << indent << "<eventsRec>" << std::endl;
  bool previousWasTransfer = previousEvent->type 
    == ReconciliationEventType::EVENT_T || previousEvent->type == ReconciliationEventType::EVENT_TL;
  if (previousWasTransfer && geneIndex == previousEvent->rightGeneIndex && event.type 
      != ReconciliationEventType::EVENT_L) {
    auto previousEventSpeciesOut = speciesTree->nodes[previousEvent->destSpeciesNode];
    os << indent << "\t<transferBack destinationSpecies=\"" << previousEventSpeciesOut->label << "\"/>" << std::endl;
  }
  
  switch(event.type) {
  case ReconciliationEventType::EVENT_None:
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

static void recursivelySaveGeneTreeRecPhyloXML(unsigned int geneIndex, 
    corax_rtree_t *speciesTree, 
    std::vector<std::vector<Scenario::Event> > &geneToEvents,
    const Scenario::Event *previousEvent,
    std::string &indent,
    ParallelOfstream &os)
{
  auto &events = geneToEvents[geneIndex];
  for (unsigned int i = 0; i < events.size() - 1; ++i) {
    os << indent << "<clade>" << std::endl;
    indent += "\t";
    auto &event = events[i];
    os << indent << "<name>" << "NULL" << "</name>" << std::endl;
    writeEventRecPhyloXML(geneIndex, speciesTree, event, previousEvent, indent, os);  
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
      writeEventRecPhyloXML(geneIndex, speciesTree, loss, previousEvent, indent, os);
      indent.pop_back();
      os << indent << "</clade>" << std::endl;
    } else {
      assert(false); 
    }
  }

  os << indent << "<clade>" << std::endl;
  indent += "\t";
  Scenario::Event &event = geneToEvents[geneIndex].back();
  auto label = event.label.size() ? event.label : "NULL";
  os << indent << "<name>" << label << "</name>" << std::endl;
  writeEventRecPhyloXML(geneIndex, speciesTree, event, previousEvent, indent, os);  


  if (!event.isLeaf()) {
    recursivelySaveGeneTreeRecPhyloXML(event.leftGeneIndex, speciesTree, geneToEvents, &event, indent, os);
    recursivelySaveGeneTreeRecPhyloXML(event.rightGeneIndex, speciesTree, geneToEvents, &event, indent, os);
  }
  for (unsigned int i = 0; i < events.size() - 1; ++i) {
    indent.pop_back();
    os << indent << "</clade>" << std::endl;
  }
  indent.pop_back();
  os << indent << "</clade>" << std::endl;
}

static void saveGeneTreeRecPhyloXML(unsigned int geneIndex,
    corax_rtree_t *speciesTree,
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  os << "<recGeneTree>" << std::endl;
  os << "<phylogeny rooted=\"true\">" << std::endl;
  std::string indent;
  Scenario::Event previousEvent;
  previousEvent.type = ReconciliationEventType::EVENT_None;
  recursivelySaveGeneTreeRecPhyloXML(geneIndex, speciesTree, geneToEvents, &previousEvent, indent, os); 
  os << "</phylogeny>" << std::endl;
  os << "</recGeneTree>" << std::endl;
}

void ReconciliationWriter::saveReconciliationRecPhyloXML(corax_rtree_t *speciesTree, 
    unsigned int geneIndex, 
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  os << "<recPhylo " << std::endl;
  os << "\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << std::endl;
  os << "\txsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\"" << std::endl;
  os << "\txmlns=\"http://www.recg.org\">" << std::endl;
  saveSpeciesTreeRecPhyloXML(speciesTree, os);
  saveGeneTreeRecPhyloXML(geneIndex, speciesTree, geneToEvents, os);
  os << "</recPhylo>";

}

static void recursivelySaveReconciliationsNewickEvents(
    corax_unode_t *node, 
    unsigned int depth, 
    std::vector<std::vector<Scenario::Event> > &geneToEvents, 
    ParallelOfstream &os)
{
  
  if(node->next) {
    auto left = node->next->back;
    auto right = node->next->next->back;
    if (depth == 0) {
      left = node->next;
      right = node->next->back;
    }
    os << "(";
    recursivelySaveReconciliationsNewickEvents(
        left, 
        depth + 1, 
        geneToEvents, 
        os);
    os << ",";
    recursivelySaveReconciliationsNewickEvents(
        right, 
        depth + 1, 
        geneToEvents, 
        os);
    os << ")";
  } 
  if (!node->next) {
    os << node->label;
  } else {
    os << Enums::getEventName(geneToEvents[node->node_index].back().type); 
  }
  if (depth == 1) {
    // we split the length of the virtual root in two
    // for each branch under the root
    os << ":" << node->length / 2.0;
  }
  if (depth > 1) {
    os << ":" << node->length;
  }
}
  
void ReconciliationWriter::saveReconciliationNewickEvents(corax_unode_t *geneRoot, 
      unsigned int virtualRootIndex,
      std::vector<std::vector<Scenario::Event> > &geneToEvent, 
      ParallelOfstream &os)
{
  corax_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = virtualRootIndex;
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  recursivelySaveReconciliationsNewickEvents(&virtualRoot, 
      0, 
      geneToEvent, 
      os);
  os << ";";

}
