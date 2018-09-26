//
// File GlobalGraph.h
// Created by: Thomas Bigot
//             Laurent Gueguen
// Last modification : vendredi 4 novembre 2016, à 10h 19
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide utilitary
   classes. This file belongs to the Bio++ Project.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _GLOBAL_GRAPH_H_
#define _GLOBAL_GRAPH_H_

#include "../Clonable.h"

#include <set>
#include <map>
#include <string>
#include "Graph.h"

namespace bpp
{
class GlobalGraph :
  public virtual Graph,
  public virtual Clonable
{
public:
  typedef Graph::NodeId Node;
  typedef Graph::EdgeId Edge;

  /**
   * The node structure type
   * Node -> ("toNodes" [DestNode,Edge],"fromNodes" [DestNode,Edge])
   * directed example: (N1)-E1->(N2)-E2->(N3) is coded as
   *     N1 -> ((N2:E1),())
   *     N2 -> ((N3:E3),(N1:E1))
   *     N3 -> ((),(N2:E2))
   * undirected example: (N1)-E1-(N2)-E2->(N3) is coded as
   *     N1 -> ((N2:E1),(N2:E1))
   *     N2 -> ((N1:E1, N3:E3),(N1:E1, N3:E3))
   *     N3 -> ((N2:E2),(N2:E2))
   */
  typedef std::map<Node, std::pair<std::map<Node, Edge>, std::map<Node, Edge> > > nodeStructureType;

  /**
   * The edge structure type
   * directed example: N1--E1-->N2 is coded as E1 -> (N1,N2)
   * undirected example: N1--E1--N2 is coded as E1 -> (N1,N2)
   */
  typedef std::map<Edge, std::pair<Node, Node> > edgeStructureType;

private:
  /**
   * is the graph directed
   */
  bool directed_;

  /**
   * List of all the subscribers.
   */
  std::set<GraphObserver*> observers_;

  /**
   * Highest used available ID for a Node.
   */
  Node highestNodeID_;
  /**
   * Highest used available ID for an Edge.
   */
  Edge highestEdgeID_;

  /**
   * Nodes and their relations.
   * see nodeStructureType documentation
   */

  nodeStructureType nodeStructure_;

  /**
   * Edges and their relations in the forward direction..
   * see edgeStructureType documentation
   */
  edgeStructureType edgeStructure_;

  /**
   * Usualy the first node of a graph. Used for algorithmic purposes.
   */
  Node root_;

  /**
   * Some types of Graphs need to know if they have been modified
   * But for a Graph, it does nothing.
   */
  virtual void topologyHasChanged_() const
  {
    // do nothing: a Graph does not care to be modified
  }

  /**
   * Tell all the observers to get the last updates.
   * Calls the method update of all the subscribers.
   */
  void notify_();

  /**
   * Creates a link between two existing nodes. If directed graph: nodeA -> nodeB.
   * Private version of link, does not check for the reciprocity.
   * Mainly called by link().
   * @param nodeA source node
   * @param nodeB target node
   * @param edge the ID of the relation
   */
  void linkInNodeStructure_(const Node& nodeA, const Node& nodeB, const Edge& edge);

  /**
   * Creates a link between two existing nodes in the edge structure.
   * If directed graph: nodeA -> nodeB.
   * Mainly called by link().
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   * @param edge the ID of the relation
   */
  void linkInEdgeStructure_(const Node& nodeA, const Node& nodeB, const Edge& edge);


  /**
   * Erase a link between two existing nodes. If directed graph: nodeA -> nodeB.
   * Private version of unLink, does not check for the reciprocity.
   * Mainly called by unLink().
   * @param nodeA source node
   * @param nodeB target node
   * @return the ID of the erased relation
   */

  Edge unlinkInNodeStructure_(const Node& nodeA, const Node& nodeB);

  /**
   * Erase a link between two existing nodes in the Edge structure.
   * Mainly called by unLink().
   * @param edge the edge to unregister
   */

  void unlinkInEdgeStructure_(const Edge& edge);

protected:
  /**
   * get the Highest Node ID (for vector sizing)
   */
  Node getHighestNodeID() const;

  /**
   * get the Highest Node ID (for vector sizing)
   */
  Edge getHighestEdgeID() const;


  /**
   * Check that a node exists. If not, throw an exception.
   * @param node node that has to be checked
   * @param name common name to give to the user in case of failure (eg: "first node")
   */
  void nodeMustExist_(const Node& node, std::string name = "") const;

  /**
   * Check that a edge exists. If not, throw an exception.
   * @param edge edge that has to be checked
   * @param name common name to give to the user in case of failure (eg: "first node")
   */
  void edgeMustExist_(const Edge& edge, std::string name = "") const;

private:
  /**
   * Private version of getIncomingNeighbors or getOutgoingNeighbors.
   * Common code of these function shared here.
   * @param node node to  in or outgoing neighbors
   * @param outgoing boolean: if true, outgoing; else incoming
   */
  std::vector<Node> getNeighbors_(const Node& node, bool outgoing = true) const;

  /**
   * Private version of getIncomingEdges or getOutgoingEdges.
   * Common code of these function shared here.
   * @param node node to  in or outgoing edges
   * @param outgoing boolean: if true, outgoing; else incoming
   */
  std::vector<Edge> getEdges_(const Node& node, bool outgoing = true) const;

  /**
   * Separate a node from all its neighbors.
   * @param node node to isolate
   */

  void isolate_(Node& node);

  /**
   * Get leaves from a starting node, filling a vector (private version).
   * @param startingNode root node
   * @param foundLeaves a vector containing all the found leaves
   * @param originNode the node where we come from, not to explore
   * @param maxRecursions  maximum number of recursion steps
   */

  void fillListOfLeaves_(const Node& startingNode, std::vector<Node>& foundLeaves, const Node& originNode, unsigned int maxRecursions) const;

  /**
   * Check that nodes are only met once to define if the graph is cyclic.
   * @param node the node to explore
   * @param metNodes a set containing all the nodes we met
   * @param originNode the node where we come from, not to explore
   */
  bool nodesAreMetOnlyOnce_(const Node& node, std::set<Node>& metNodes, const Node& originNode) const;

  /**
   * output a node to DOT format (recursive)
   */

  void nodeToDot_(const Node& node, std::ostream& out, std::set<std::pair<Node, Node> >& alreadyFigured) const;

public:
  /** @name General Management
   *  Misc & constructors
   */
  // /@{


  /**
   * Constructor
   * @param directed true if the graph is directed.
   */
  GlobalGraph(bool directed = false);

  GlobalGraph(const GlobalGraph& gg);

  GlobalGraph& operator=(const GlobalGraph& gg);

  GlobalGraph* clone() const {return new GlobalGraph(*this); }

  ~GlobalGraph() {}

protected:
  /**
   * set the root node to an existing node. Will not affect the topology.
   * @param newRoot the new root
   */

  void setRoot(Graph::NodeId newRoot);

public:
  /**
   * get the root node
   */

  Graph::NodeId getRoot() const;

  /**
   * Make the graph directed
   * - changes the property
   * - de-duplicate the relations:
   *  eg: A - B in undirected is represented as A->B and B->A
   *  in directed, becomes A->B only
   *
   * Please note that the resulting directions are totaly arbritrary.
   * One might consider to use the makeLocalRoot method.
   */
  void makeDirected();

  /**
   * Make the graph directed
   * - changes the property
   * - de-duplicate the relations:
   *    eg: A - B in directed is represented as A->B
   *        in undirected, becomes A->B and B->A
   * If the directed graph already contains reciprocal relations,
   * such as A->B and B->A, the method will throw an exception.
   */
  void makeUndirected();

  // /@}


  /** @name Relations management
   *  Modificating the structure of the graph.
   */
  // /@{

  /**
   * Creates an orphaned node.
   * @return the new node
   */

  Graph::NodeId createNode();

  /**
   * Creates a node linked to an existing node.
   * @param origin existing node. In a directed graph: origin -> newNode.
   * @return the new node
   */

  Graph::NodeId createNodeFromNode(Graph::NodeId origin);

  /**
   * Creates new node on an existing Edge. A -> B will be A -> N -> B
   * @param edge existing edge.
   * @return the new node
   */

  Graph::NodeId createNodeOnEdge(Graph::EdgeId edge);

  /**
   * Creates a node linked to new node, splitting an edge.
   * @param origin existing edge. In a directed graph: origin -> newNode.
   * @return the new node
   */

  Graph::NodeId createNodeFromEdge(Graph::NodeId origin);

protected:
  /**
   * Creates a link between two existing nodes. If directed graph: nodeA -> nodeB.
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   * @return the new edge
   */

  Graph::EdgeId link(Graph::NodeId nodeA, Graph::NodeId nodeB);

  /**
   * Sets a link between two existing nodes, using existing edge. If
   * directed graph: nodeA -> nodeB.
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   * @param edgeID the used edge
   */

  void link(Graph::NodeId nodeA, Graph::NodeId nodeB, GlobalGraph::Edge edgeID);

  /**
   * Switch the edge  between two existing nodes.
   *
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   */

  void switchNodes(Graph::NodeId nodeA, Graph::NodeId nodeB);


  /**
   * Remove all links between two existing nodes. If directed graph: nodeA -> nodeB.
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   * @return vector of deleted edges
   */

  std::vector<Graph::EdgeId> unlink(Graph::NodeId nodeA, Graph::NodeId nodeB);

public:
  /**
   * Delete one node
   * @param node node to be deleted
   */

  void deleteNode(Graph::NodeId node);

  // /@}


  /** @name Observers Management
   *  Managing communication with the observers: subscribe, unsubscribe.
   */
  // /@{

  /**
   * Attach a new observer to this Graph.
   * As a subscriber, the observer will be warned of all the changes.
   */
  void registerObserver(GraphObserver* observer);
  /**
   * Detach an observer from this Graph.
   * The observer will not be warned of changes anymore.
   */
  void unregisterObserver(GraphObserver* observer);
  // /@}


  /** @name Nodes Functions
   *  These methodes of the graph concern the node management.
   */
  // /@{

  template<typename T, bool is_const>
  friend class NodesIteratorClass;

  template<typename T, bool is_const>
  friend class EdgesIteratorClass;

  /*
   * @brief builds iterator on all Nodes of the graph
   *
   */

  std::unique_ptr<Graph::NodeIterator> allNodesIterator();
  std::unique_ptr<Graph::NodeIterator> allNodesIterator() const;


  /*
   * @brief builds iterator on the neighbor nodes of a Node
   *
   */

  std::unique_ptr<Graph::NodeIterator> outgoingNeighborNodesIterator(NodeId node);
  std::unique_ptr<Graph::NodeIterator> outgoingNeighborNodesIterator(NodeId node) const;

  /*
   * @brief builds iterator on the neighbor nodes of a Node
   *
   */

  std::unique_ptr<Graph::NodeIterator> incomingNeighborNodesIterator(NodeId node);
  std::unique_ptr<Graph::NodeIterator> incomingNeighborNodesIterator(NodeId node) const;

  /**
   * Get the number of nodes in the graph.
   */

  size_t getNumberOfNodes() const;

  /**
   * Get the number of edges in the graph.
   */

  size_t getNumberOfEdges() const;

  /**
   * Get the degree of a node (ie the number of neighbors) in the graph.
   * @param node the node one wants to count its neighbors
   * @return the number of neighbors
   */

  size_t getDegree(Graph::NodeId node) const;

  /**
   * Says if  a node is a leaf (ie has at most one neighbor).
   */

  bool isLeaf(Graph::NodeId node) const;

  /**
   * Get the number of  neighbors  of a node in the graph.
   * @param node the node one wants to count its neighbors
   * @return the number of neighbors
   */

  size_t getNumberOfNeighbors(Graph::NodeId node) const;

  /**
   * Get the number of outgoing neighbors  of a node (ie the number of sons) in the graph.
   * @param node the node one wants to count its sons
   * @return the number of outgoing neighbors
   */

  size_t getNumberOfOutgoingNeighbors(Graph::NodeId node) const;

  /**
   * Get the number of incoming neighbors  of a node (ie the number of fathers) in the graph.
   * @param node the node one wants to count its fathers
   * @return the number of incoming neighbors
   */

  size_t getNumberOfIncomingNeighbors(Graph::NodeId node) const;

  /**
   * Get all the neighbors of a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the neighbors
   */

  std::vector<Graph::NodeId> getNeighbors(Graph::NodeId node) const;

  /**
   * In an directed graph, get all the neighbors which
   * are leaving a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the outgoing neighbors
   */
  std::vector<Graph::NodeId> getOutgoingNeighbors(Graph::NodeId node) const;

  /**
   * In an directed graph, get all the neighbors which
   * are coming to a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the incoming neighbors
   */
  std::vector<Graph::NodeId> getIncomingNeighbors(Graph::NodeId node) const;

  /**
   * Get the leaves of a graph, ie, nodes with only one neighbor,
   * starting from a peculiar node.
   * @param node the starting node
   * @param maxDepth the maximum number of allowed depth.
   * @return a vector containing the leaves
   */
  std::vector<Graph::NodeId> getLeavesFromNode(Graph::NodeId node, unsigned int maxDepth) const;

  /**
   * Get all leaves of a graph, ie, nodes with no son (or only one
   * neighbor is not directet).
   * @return a vector containing the leaves
   */

  std::vector<Graph::NodeId> getAllLeaves() const;
  std::set<NodeId> getSetOfAllLeaves() const;

  /**
   * Get all the nodes.
   * @return a vector containing the  nodes
   */

  std::vector<Graph::NodeId> getAllNodes() const;

  /**
   * Get all the inner nodes, ie, nodes with degree > 1.
   * @return a vector containing the inner nodes
   */

  std::vector<Graph::NodeId> getAllInnerNodes() const;

  /**
   * Get nodes located at the extremities of an edge
   *
   * @return a pair of the IDs of the Nodes at each extremity of the edge
   *        example : N1--E1-->N2; getNodes(E1) will return (N1,N2);
   */
  std::pair<Graph::NodeId, Graph::NodeId> getNodes(Graph::EdgeId edge) const;

  /**
   * Get node located at the top of an edge
   *
   * @return  the Node at the top the edge
   *        example : N1--E1-->N2; getTop(E1) will return N1;
   */

  Graph::NodeId getTop(Graph::EdgeId edge) const;

  /**
   * Get node located at the bottom of an edge
   *
   * @return  the Node at the bottom the edge
   *        example : N1--E1-->N2; getBottom(E1) will return N2;
   */

  Graph::NodeId getBottom(Graph::EdgeId edge) const;


  // /@}

  /** @name Topological Properties
   *  These methodes check some topological properties.
   */
  // /@{

  /**
   * Is the graph a tree?
   * @return false if a node is met more than one time browsing the graph
   */

  bool isTree() const;

  /**
   * Is the graph directed acyclic?
   * @return true if an edge is met more than one time browsing the graph
   */

  bool isDA() const;

  /**
   * Orientates the graph hanging from the root
   */

  void orientate();

  /**
   * Is the graph directed?
   * @return true the type of the graph is directed
   */
  bool isDirected() const;


  /**
   * Does the graph contain reciprocal relations such as A->B and B->A?
   * @return true if one of them is seen in the structure
   */
  bool containsReciprocalRelations() const;


  // /@}

  /*
   * @brief builds iterator on all Nodes of the graph
   *
   */

  std::unique_ptr<EdgeIterator> allEdgesIterator();
  std::unique_ptr<EdgeIterator> allEdgesIterator() const;

  /*
   * @brief builds iterator on the outgoing neighbor nodes of a Node
   *
   */

  std::unique_ptr<EdgeIterator> outgoingEdgesIterator(NodeId node);
  std::unique_ptr<EdgeIterator> outgoingEdgesIterator(NodeId node) const;

  /*
   * @brief builds iterator on the incoming neighbor nodes of a Node
   *
   */

  std::unique_ptr<EdgeIterator> incomingEdgesIterator(NodeId node);
  std::unique_ptr<EdgeIterator> incomingEdgesIterator(NodeId node) const;

  /**
   * Get all the edges to/from a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the edges
   */

  std::vector<Graph::EdgeId> getEdges(Graph::NodeId node) const;

  /**
   * In an directed graph, get all the edges which
   * are leaving a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the outgoing edges
   */
  std::vector<Graph::EdgeId> getOutgoingEdges(Graph::NodeId node) const;

  /**
   * In an directed graph, get all the edges which
   * are coming to a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the incoming edges
   */
  std::vector<Graph::EdgeId> getIncomingEdges(Graph::NodeId node) const;

  /**
   * Returns the Edge between two nodes
   * @param nodeA if directed, origin node
   * @param nodeB if directed, destination node
   * @return the edge between these two nodes
   */

  Graph::EdgeId getEdge(Graph::NodeId nodeA, Graph::NodeId nodeB) const;

  /**
   * Returns the Edge between two nodes, trying both directions
   * @param nodeA any node implied in the relation
   * @param nodeB any other node implied in the relation
   * @return the edge between these two nodes
   */

  Graph::EdgeId getAnyEdge(Graph::NodeId nodeA, Graph::NodeId nodeB) const;

  /**
   * Get all edges of a graph.
   * @return a vector containing the edges
   */

  std::vector<Graph::EdgeId> getAllEdges() const;

  // /@}


  /** @name Updating the changes on the observers
   *  These methodes aim to trigger some changes to the observers
   */
  // /@{

  /**
   * Trigger E objects deleting on the observers
   * @param edgesToDelete list of edges to delete
   */

  void notifyDeletedEdges(const std::vector<Graph::EdgeId>& edgesToDelete) const;

  /**
   * Trigger N objects deleting on the observers
   * @param nodesToDelete list of edges to delete
   */
  void notifyDeletedNodes(const std::vector<Graph::NodeId>& nodesToDelete) const;


  // /@}

  /**
   * Output the graph in DOT format
   * @param out a ostream where the DOT format will be output
   * @param name a string naming the graph
   */

  void outputToDot(std::ostream& out, const std::string& name) const;

  template<class N, class E, class GraphImpl>
  friend class AssociationGraphImplObserver;
};


/************************************************/
/************************************************/
/************************************************/
/* ITERATORS */
/************************************************/

/************************************************/
/* NODES ITERATORS */
/************************************************/

template<class T, bool is_const>
class NodesIteratorClass :
  virtual public Graph::NodeIterator
{
  ~NodesIteratorClass<T, is_const>(){}
};


template<bool is_const>
class NodesIteratorClass<Graph::ALLGRAPHITER, is_const> :
  virtual public Graph::NodeIterator
{
private:
  typedef typename std::conditional<is_const,
                                    GlobalGraph::nodeStructureType::const_iterator,
                                    GlobalGraph::nodeStructureType::iterator >::type itType;

  itType it_,  begin_, end_;

public:
  template<bool B = is_const>
  NodesIteratorClass<Graph::ALLGRAPHITER, is_const>(const GlobalGraph &gg, typename std::enable_if<B>::type * = 0) : it_(gg.nodeStructure_.begin()),
    begin_(gg.nodeStructure_.begin()),
    end_(gg.nodeStructure_.end()) {}

  template<bool B = is_const>
  NodesIteratorClass<Graph::ALLGRAPHITER, is_const>(GlobalGraph & gg, typename std::enable_if<!B>::type * = 0) : it_(gg.nodeStructure_.begin()),
    begin_(gg.nodeStructure_.begin()),
    end_(gg.nodeStructure_.end()) {}

  ~NodesIteratorClass<Graph::ALLGRAPHITER, is_const>(){}

  void next() { it_++; }
  bool end() const { return it_ == end_;  }
  void start() { it_ = begin_; }

  Graph::NodeId operator*() {return it_->first; }
};


/**
 * @brief Abstract class for neighbor iterators
 *
 **/

template<bool is_const>
class NeighborIteratorClass
{
protected:
  typedef typename std::conditional<is_const,
                                    const std::map<GlobalGraph::Node, GlobalGraph::Edge>&,
                                    std::map<GlobalGraph::Node, GlobalGraph::Edge>&>::type mapType;

  mapType map_;

  typedef typename std::conditional<is_const,
                                    std::map<GlobalGraph::Node, GlobalGraph::Edge>::const_iterator,
                                    std::map<GlobalGraph::Node, GlobalGraph::Edge>::iterator>::type itType;

  itType it_, begin_, end_;

public:
  virtual ~NeighborIteratorClass<is_const>(){}

  template<bool B = is_const>
  NeighborIteratorClass<is_const>(const std::map<GlobalGraph::Node, GlobalGraph::Edge> &map, typename std::enable_if<B>::type * = 0) :
    map_(map),
    it_(map_.begin()),
    begin_(map_.begin()),
    end_(map_.end()) {}

  template<bool B = is_const>
  NeighborIteratorClass<is_const>(std::map<GlobalGraph::Node, GlobalGraph::Edge> &map, typename std::enable_if<!B>::type * = 0) :
    map_(map),
    it_(map_.begin()),
    begin_(map_.begin()),
    end_(map_.end()) {}

  void next() { it_++; }
  bool end() const { return it_ == end_;  }
  void start() { it_ = begin_; }
};


/**
 * @brief
 *
 **/

template<bool is_const>
class NodesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const> :
  virtual public NeighborIteratorClass<is_const>,
  virtual public Graph::NodeIterator
{
public:
  NodesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const>(const GlobalGraph &gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.first) {}

  NodesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const>(GlobalGraph & gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.first) {}

  ~NodesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const>(){}

  void next() { NeighborIteratorClass<is_const>::next(); }
  bool end() const { return NeighborIteratorClass<is_const>::end(); }
  void start() { NeighborIteratorClass<is_const>::start(); }

  Graph::NodeId operator*() {return NeighborIteratorClass<is_const>::it_->first; }
};


template<bool is_const>
class NodesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const> :
  virtual public NeighborIteratorClass<is_const>,
  virtual public Graph::NodeIterator
{
public:
  NodesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const>(const GlobalGraph &gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.second) {}

  NodesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const>(GlobalGraph & gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.second) {}

  ~NodesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const>(){}

  void next() { NeighborIteratorClass<is_const>::next(); }
  bool end() const { return NeighborIteratorClass<is_const>::end(); }
  void start() { NeighborIteratorClass<is_const>::start(); }

  Graph::NodeId operator*() {return NeighborIteratorClass<is_const>::it_->first; }
};


/************************************************/
/* EDGES ITERATORS */
/************************************************/

template<class T, bool is_const>
class EdgesIteratorClass :
  virtual public Graph::EdgeIterator
{};

template<bool is_const>
class EdgesIteratorClass<Graph::ALLGRAPHITER, is_const> :
  virtual public Graph::EdgeIterator
{
private:
  typedef typename std::conditional<is_const,
                                    GlobalGraph::edgeStructureType::const_iterator,
                                    GlobalGraph::edgeStructureType::iterator >::type itType;

  itType it_, begin_, end_;

public:
  template<bool B = is_const>
  EdgesIteratorClass<Graph::ALLGRAPHITER, is_const>(const GlobalGraph &gg, typename std::enable_if<B>::type * = 0) : it_(gg.edgeStructure_.begin()),
    begin_(gg.edgeStructure_.begin()),
    end_(gg.edgeStructure_.end()) {}

  template<bool B = is_const>
  EdgesIteratorClass<Graph::ALLGRAPHITER, is_const>(GlobalGraph & gg, typename std::enable_if<!B>::type * = 0) : it_(gg.edgeStructure_.begin()),
    begin_(gg.edgeStructure_.begin()),
    end_(gg.edgeStructure_.end()) {}

  ~EdgesIteratorClass<Graph::ALLGRAPHITER, is_const>(){}

  void next() { it_++; }
  bool end() const { return it_ == end_;  }
  void start() { it_ = begin_; }

  Graph::EdgeId operator*() {return it_->first; }
};


template<bool is_const>
class EdgesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const> :
  public NeighborIteratorClass<is_const>,
  public Graph::EdgeIterator
{
public:
  EdgesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const>(const GlobalGraph &gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.first) {}

  EdgesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const>(GlobalGraph & gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.first) {}

  ~EdgesIteratorClass<Graph::OUTGOINGNEIGHBORITER, is_const>(){}

  void next() { NeighborIteratorClass<is_const>::next(); }
  bool end() const { return NeighborIteratorClass<is_const>::end(); }
  void start() { NeighborIteratorClass<is_const>::start(); }

  Graph::EdgeId operator*() {return NeighborIteratorClass<is_const>::it_->second; }
};

template<bool is_const>
class EdgesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const> :
  public NeighborIteratorClass<is_const>,
  public Graph::EdgeIterator
{
public:
  EdgesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const>(const GlobalGraph &gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.second) {}

  EdgesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const>(GlobalGraph & gg, GlobalGraph::NodeId node) : NeighborIteratorClass<is_const>(gg.nodeStructure_.find(node)->second.second) {}

  ~EdgesIteratorClass<Graph::INCOMINGNEIGHBORITER, is_const>(){}

  void next() { NeighborIteratorClass<is_const>::next(); }
  bool end() const { return NeighborIteratorClass<is_const>::end(); }
  void start() { NeighborIteratorClass<is_const>::start(); }

  Graph::EdgeId operator*() {return NeighborIteratorClass<is_const>::it_->second; }
};
}

#endif
