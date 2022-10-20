#pragma once

#include <vector>
#include <string>
#include <util/enums.hpp>
#include <util/types.hpp>
#include <memory>
#include <unordered_set>
#include <corax/corax.h>
#include <trees/PLLRootedTree.hpp>

/*
typedef struct corax_utree_s corax_utree_t;
typedef struct corax_unode_s corax_unode_t;
typedef struct corax_rtree_s corax_rtree_t;
typedef struct corax_rnode_s corax_rnode_t;
*/
class ParallelOfstream;
typedef std::unordered_set<std::string> OrthoGroup;
typedef std::shared_ptr<OrthoGroup> OrthoGroupPtr;
typedef std::vector<OrthoGroupPtr> OrthoGroups;

struct SpeciesEvents {
  unsigned int LeafCount;
  unsigned int DCount;
  unsigned int SCount;
  unsigned int SLCount;
  unsigned int TCount;
  unsigned int TLCount;
  SpeciesEvents():
    LeafCount(0), 
    DCount(0), 
    SCount(0),
    SLCount(0),
    TCount(0),
    TLCount(0)
  {}


  double speciesFrequency() const {
    return SCount + TCount + DCount + SLCount + TLCount;
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
      destSpeciesNode(INVALID_NODE_ID),
      pllTransferedGeneNode(nullptr),
      pllDestSpeciesNode(nullptr)
    {}
    ReconciliationEventType type;
    unsigned int geneNode;
    unsigned int speciesNode;
    unsigned int destSpeciesNode;     // for transfers only
   
    // left and gene indices:
    // - for speciation: left gene goes to left species
    // - for transfer: left gene goes to source species and right gene 
    // goes to recieving species
    // - for duplication, the order doesn't matter
    // - for other events, those members are irrelevant
    unsigned int leftGeneIndex;
    unsigned int rightGeneIndex;

    // temporary variables for event inference
    corax_unode_t *pllTransferedGeneNode;
    corax_rnode_t *pllDestSpeciesNode;
    
    std::string label;

    /**
     *  Is this event valid? (if not, something went wrong)
     */
    bool isValid() const { return speciesNode != INVALID_NODE_ID; }
    bool isLeaf() const { return type == ReconciliationEventType::EVENT_None;}
  };


  /**
   * Default constructor
   */
  Scenario(): 
    _eventsCount(static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid), 0), 
    _geneRoot(nullptr), 
    _rootedTree(nullptr),
    _speciesTree(nullptr),
    _virtualRootIndex(INVALID_NODE_ID) 
  {}

  // forbid copy
  Scenario(const Scenario &) = delete;
  Scenario & operator = (const Scenario &) = delete;
  Scenario(Scenario &&) = delete;
  Scenario & operator = (Scenario &&) = delete;
 
  void setGeneRoot(corax_unode_t *geneRoot) {_geneRoot = geneRoot;}
  void setSpeciesTree(PLLRootedTree *speciesTree) {
    _rootedTree = speciesTree;
    _speciesTree = _rootedTree->getRawPtr();
  }
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
    unsigned int destSpeciesNode);

  /**
   * Various methods to save information from the Scenario
   */
  void saveEventsCounts(const std::string &filename, bool masterRankOnly = true); 
  void saveTransfers(const std::string &filename, bool masterRankOnly = true); 
  void saveReconciliation(const std::string &filename, ReconciliationFormat format, bool masterRankOnly = true);
  void saveReconciliation(ParallelOfstream &os, ReconciliationFormat format);
  static void saveTransferPairCountGlobal(PLLRootedTree &speciesTree,
      std::vector<Scenario> &scenarios,
      const std::string &filename);

  void saveLargestOrthoGroup(std::string &filename, bool masterRankOnly = true) const;
  void saveAllOrthoGroups(std::string &filename, bool masterRankOnly = true) const;
  void savePerSpeciesEventsCounts(const std::string &filename, bool masterRankOnl = true);
  void gatherReconciliationStatistics(PerSpeciesEvents &perSpeciesEvents) const;
  void countTransfers(const StringToUint &labelToId,
      MatrixUint &count);

  unsigned int getEventCount(ReconciliationEventType type) const {
    return _eventsCount[static_cast<unsigned int>(type)];
  }
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


  std::vector<corax_unode_t *> &getGeneNodeBuffer() {return _geneNodeBuffer;}
  corax_unode_t *generateVirtualGeneRoot();
  void generateGeneChildren(corax_unode_t *geneNode, 
      corax_unode_t *&leftGeneNode,
      corax_unode_t *&rightGeneNode);


  corax_unode_t *getGeneRoot() const {return _geneRoot;}
  unsigned int getVirtualRootIndex() const { return _virtualRootIndex;}
  corax_rtree_t *getSpeciesTree() const {return _speciesTree;}
  const std::vector<std::vector<Event> > &
    getGeneIdToEvents() const {return _geneIdToEvents;}
private:
  static const char *eventNames[];
  std::vector<Event> _events;
  std::vector<unsigned int> _eventsCount;
  std::vector<std::vector<Event> > _geneIdToEvents;
  corax_unode_t *_geneRoot;
  PLLRootedTree *_rootedTree;
  corax_rtree_s *_speciesTree;
  unsigned int _virtualRootIndex;
  std::vector<corax_unode_t *> _geneNodeBuffer;
  typedef std::vector< std::vector <int> > ScenarioBlackList;
  std::unique_ptr<ScenarioBlackList> _blacklist;
  OrthoGroup *getLargestOrthoGroupRec(corax_unode_t *geneNode, bool isVirtualRoot) const;
  void getAllOrthoGroupRec(corax_unode_t *geneNode,
      OrthoGroups &orthogroups,
      OrthoGroupPtr &currentOrthoGroup,
      bool isVirtualRoot) const;

};



