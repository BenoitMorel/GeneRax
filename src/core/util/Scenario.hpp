#pragma once

#include <vector>
#include <string>
#include <util/enums.hpp>
#include <memory>
#include <unordered_set>
extern "C" {
#include <pll.h>
}

typedef struct pll_utree_s pll_utree_t;
typedef struct pll_unode_s pll_unode_t;
typedef struct pll_rtree_s pll_rtree_t;
typedef struct pll_rnode_s pll_rnode_t;
class ParallelOfstream;
typedef std::unordered_set<std::string> OrthoGroup;
typedef std::shared_ptr<OrthoGroup> OrthoGroupPtr;
typedef std::vector<OrthoGroupPtr> OrthoGroups;

struct SpeciesEvents {
  unsigned int DCount;
  unsigned int SCount;
  unsigned int SLCount;
  unsigned int TCount;
  unsigned int TLCount;
  SpeciesEvents():
    DCount(0), 
    SCount(0),
    SLCount(0),
    TCount(0),
    TLCount(0)
  {}


  double speciesFrequency() const {
    return SCount + TCount;
  }
};

struct PerSpeciesEvents {
  std::vector<SpeciesEvents> events;
  PerSpeciesEvents(unsigned int speciesNodeCount = 0):
    events(speciesNodeCount) 
  {}
  void parallelSum();
};

/*
 * Store the set of events that reconciles a gene tree with a species tree
 */
class Scenario {
public: 
  // value representing an invalid node ID (0 is a valid node ID)
  static const unsigned int INVALID_NODE_ID = static_cast<unsigned int>(-1);
 
  /**
   *  Reconciliation event description: type of event, and
   *  involved gene and species nodes
   */
  struct Event {
    Event(): 
      type(ReconciliationEventType::EVENT_S), 
      geneNode(INVALID_NODE_ID), 
      speciesNode(INVALID_NODE_ID), 
      transferedGeneNode(INVALID_NODE_ID), 
      destSpeciesNode(INVALID_NODE_ID),
      pllTransferedGeneNode(nullptr),
      pllDestSpeciesNode(nullptr),
      cross(false)
    {}
    ReconciliationEventType type;
    unsigned int geneNode;
    unsigned int speciesNode;
    unsigned int transferedGeneNode;  // for transfers only 
    unsigned int destSpeciesNode;     // for transfers only
    
    // temporary variables for event inference
    pll_unode_t *pllTransferedGeneNode;
    pll_rnode_t *pllDestSpeciesNode;
    bool cross; // used to decide which gene (left/right) goes to which species (left/right)
                // if set to true, left goes to right and right goes to left
                // if set to false, left goes to left and right goes to right
    
   
    /**
     *  Is this event valid? (if not, something went wrong)
     */
    bool isValid() const { return speciesNode != INVALID_NODE_ID; }
  };


  /**
   * Default constructor
   */
  Scenario(): 
    _eventsCount(static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid), 0), 
    _geneRoot(nullptr), 
    _virtualRootIndex(INVALID_NODE_ID) 
  {}

  // forbid copy
  Scenario(const Scenario &) = delete;
  Scenario & operator = (const Scenario &) = delete;
  Scenario(Scenario &&) = delete;
  Scenario & operator = (Scenario &&) = delete;
 
  void setGeneRoot(pll_unode_t *geneRoot) {_geneRoot = geneRoot;}
  void setSpeciesTree(pll_rtree_t *speciesTree) {_speciesTree = speciesTree;}
  void setVirtualRootIndex(unsigned int virtualRootIndex) {_virtualRootIndex = virtualRootIndex;}

  /**
   * Various methods to add an event in the Scnenario
   */
  void addEvent(const Event &event);
  void addEvent(ReconciliationEventType type, 
      unsigned int geneNode, 
      unsigned int speciesNode, 
      unsigned int destSpeciesNode = INVALID_NODE_ID);
  void addTransfer(ReconciliationEventType type, 
    unsigned int geneNode, 
    unsigned int speciesNode, 
    unsigned int transferedGeneNode,
    unsigned int destSpeciesNode);

  /**
   * Various methods to save information from the Scenario
   */
  void saveEventsCounts(const std::string &filename, bool masterRankOnly = true); 
  void saveTransfers(const std::string &filename, bool masterRankOnly = true); 
  void saveReconciliation(const std::string &filename, ReconciliationFormat format, bool masterRankOnly = true);
  void saveReconciliation(ParallelOfstream &os, ReconciliationFormat format);

  void saveLargestOrthoGroup(std::string &filename, bool masterRankOnly = true) const;
  void saveAllOrthoGroups(std::string &filename, bool masterRankOnly = true) const;
  void savePerSpeciesEventsCounts(const std::string &filename, bool masterRankOnl = true);
  void gatherReconciliationStatistics(PerSpeciesEvents &perSpeciesEvents) const;
  /**
   * The following methods help blacklisting pairs of gene and species nodes.
   * The motivation is that both ML and sampling reconciliation inference algorithm
   * could get stuck in an infinite loop of transfers, because of some approximations
   * in the likelihood computation. 
   * The blacklist can be used by these algorithm to detect this infinite loop. 
   */
  void initBlackList(unsigned int genesNumber, unsigned int speciesNumber);
  void blackList(unsigned int geneNode, unsigned int speciesNode);
  bool isBlacklisted(unsigned int geneNode, unsigned int speciesNode);
  void resetBlackList();
private:
  static const char *eventNames[];
  std::vector<Event> _events;
  std::vector<unsigned int> _eventsCount;
  std::vector<std::vector<Event> > _geneIdToEvents;
  pll_unode_t *_geneRoot;
  pll_rtree_t *_speciesTree;
  unsigned int _virtualRootIndex;
  typedef std::vector< std::vector <int> > ScenarioBlackList;
  std::unique_ptr<ScenarioBlackList> _blacklist;
  OrthoGroup *getLargestOrthoGroupRec(pll_unode_t *geneNode, bool isVirtualRoot) const;
  void getAllOrthoGroupRec(pll_unode_t *geneNode,
      OrthoGroups &orthogroups,
      OrthoGroupPtr &currentOrthoGroup,
      bool isVirtualRoot) const;

};



