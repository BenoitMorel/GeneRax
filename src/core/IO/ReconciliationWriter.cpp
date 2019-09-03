#include "ReconciliationWriter.hpp"

#include <IO/ParallelOfstream.hpp>
extern "C" {
#include <pll.h>
}

static void recursivelySaveReconciliationsNHX(pll_rtree_t *speciesTree,  pll_unode_t *node, std::vector<Scenario::Event> &events, ParallelOfstream &os)
{
  if(node->next) {
    os << "(";
    recursivelySaveReconciliationsNHX(speciesTree, node->next->back, events, os);
    os << ",";
    recursivelySaveReconciliationsNHX(speciesTree, node->next->next->back, events, os);
    os << ")";
  } 
  if (node->label) {
    os << node->label;
  } else {
    os << "n" << node->node_index; 
  }
  os << ":" << node->length;
  auto event = events[node->node_index];
  if (event.isValid()) {
    os << "[&&NHX";
    if (speciesTree->nodes[event.speciesNode]->label) {
      os << ":S=" << speciesTree->nodes[event.speciesNode]->label;
    }
    os << ":D=" << (event.type == EVENT_D ? "Y" : "N" );
    os << ":H=" << (event.type == EVENT_T || event.type == EVENT_TL ? "Y" : "N" );
    if (event.type == EVENT_T || event.type == EVENT_TL) {
      assert(speciesTree->nodes[event.speciesNode]->label);
      assert(speciesTree->nodes[event.destSpeciesNode]->label);
      os << "@" << speciesTree->nodes[event.speciesNode]->label;
      os << "@" << speciesTree->nodes[event.destSpeciesNode]->label;
    }
    os << ":B=" << node->length;
    os << "]";
  }
}
  
void ReconciliationWriter::saveReconciliationNHX(pll_rtree_t *speciesTree, 
    pll_unode_t *geneRoot, 
    std::vector<Scenario::Event> &events, 
    const std::string &filename, 
    bool masterRankOnly) 
{
  ParallelOfstream os(filename, masterRankOnly);
  os << "(";
  recursivelySaveReconciliationsNHX(speciesTree, geneRoot, events, os);
  os << ",";
  recursivelySaveReconciliationsNHX(speciesTree, geneRoot->back, events, os);
  os << ");";
}


static void recursivelySaveSpeciesTreeRecPhyloXML(pll_rnode_t *node, std::string &indent, ParallelOfstream &os)
{
  if (!node) {
    return;
  }
  indent += "\t";
  os << indent << "<clade>" << std::endl;
  os << indent << "\t<name>" << node->label << "</name>" << std::endl;
  recursivelySaveSpeciesTreeRecPhyloXML(node->left, indent, os);
  recursivelySaveSpeciesTreeRecPhyloXML(node->right, indent, os);
  os << indent << "</clade>" << std::endl;
  indent.pop_back();
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
    Scenario::Event &event,
    const Scenario::Event &lastEvent,
    std::string &indent, 
    ParallelOfstream &os)
{
  auto species = speciesTree->nodes[event.speciesNode];
  pll_rnode_t *speciesOut = 0;
  os << indent << "<eventsRec>" << std::endl;
  if (lastEvent.type == EVENT_T && geneTree->node_index == lastEvent.transferedGeneNode) {
    auto lastEventSpeciesOut = speciesTree->nodes[lastEvent.destSpeciesNode];
    os << indent << "\t<transferBack destinationSpecies=\"" << lastEventSpeciesOut->label << "\"/>" << std::endl;
  }
  
  switch(event.type) {
  case EVENT_None:
    assert(geneTree->next == 0);
    assert(species->left == 0 && species->right == 0);
    os << indent << "\t<leaf speciesLocation=\"" << species->label << "\"/>" <<  std::endl;
    break;
  case EVENT_S:
    os << indent << "\t<speciation speciesLocation=\"" << species->label << "\"/>" << std::endl;
    break;
  case EVENT_D:
    os << indent << "\t<duplication speciesLocation=\"" << species->label << "\"/>" << std::endl;
    break;
  case EVENT_T:
    speciesOut = speciesTree->nodes[event.speciesNode];
    os << indent << "\t<branchingOut speciesLocation=\"" << speciesOut->label << "\"/>" << std::endl;
    break; 
  default:
    const char *eventNames[]  = {"S", "SL", "D", "T", "TL", "None", "Invalid"};
    std::cerr << "please handle " << eventNames[(unsigned int)event.type] << std::endl; 
    break;
  }
  os << indent << "</eventsRec>" << std::endl;
}

static void recursivelySaveGeneTreeRecPhyloXML(pll_unode_t *geneTree, 
    pll_rtree_t *speciesTree, 
    std::vector<Scenario::Event> &events,
    const Scenario::Event &lastEvent,
    std::string &indent,
    ParallelOfstream &os)
{
  if (!geneTree) {
    return;
  }
  indent += "\t";
  os << indent << "<clade>" << std::endl;
  Scenario::Event &event = events[geneTree->node_index];
  os << indent << "<name>" << (geneTree->label ? geneTree->label : "NULL") << "</name>" << std::endl;
  writeEventRecPhyloXML(geneTree, speciesTree, event, lastEvent, indent, os);  

  if (geneTree->next) {
    recursivelySaveGeneTreeRecPhyloXML(geneTree->next->back, speciesTree, events, event, indent, os);
    recursivelySaveGeneTreeRecPhyloXML(geneTree->next->next->back, speciesTree, events, event, indent, os);
  }
  os << indent << "</clade>" << std::endl;
  indent.pop_back();
}

static void saveGeneTreeRecPhyloXML(pll_unode_t *geneTree,
    pll_rtree_t *speciesTree,
    std::vector<Scenario::Event> &events, 
    ParallelOfstream &os)
{
  os << "<recGeneTree>" << std::endl;
  os << "<phylogeny rooted=\"true\">" << std::endl;
  std::string indent;
  Scenario::Event noEvent;
  noEvent.type = EVENT_None;
  recursivelySaveGeneTreeRecPhyloXML(geneTree, speciesTree, events, noEvent, indent, os); 
  os << "</phylogeny>" << std::endl;
  os << "</recGeneTree>" << std::endl;
}

void ReconciliationWriter::saveReconciliationRecPhyloXML(pll_rtree_t *speciesTree, 
    pll_unode_t *geneRoot, 
    std::vector<Scenario::Event> &events, 
    const std::string &filename, 
    bool masterRankOnly) 
{
  ParallelOfstream os(filename, masterRankOnly);
  os << "<recPhylo " << std::endl;
  os << "\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << std::endl;
  os << "\txsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\"" << std::endl;
  os << "\txmlns=\"http://www.recg.org\">" << std::endl;
  saveSpeciesTreeRecPhyloXML(speciesTree, os);
  saveGeneTreeRecPhyloXML(geneRoot, speciesTree, events, os);
  os << "</recPhylo>";
}

