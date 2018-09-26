//
// File Graph.h
// Created by: Thomas Bigot
// Last modification : vendredi 4 novembre 2016, à 10h 21
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

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <set>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <memory>


// forward declaration to avoid circular dependancies.
// since we do not need its size in this header file (only using pointers to it)
namespace bpp
{
class GraphObserver;
}

namespace bpp
{
class Graph
{
public:
  typedef unsigned int NodeId;
  typedef unsigned int EdgeId;


  virtual ~Graph(){}

protected:
  /**
   * set the root node to an existing node. Will not affect the topology.
   * @param newRoot the new root
   */

  virtual void setRoot(NodeId newRoot) = 0;

public:
  /**
   * get the root node
   */

  virtual NodeId getRoot() const = 0;

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

  virtual void makeDirected() = 0;

  /**
   * Make the graph directed
   * - changes the property
   * - de-duplicate the relations:
   *    eg: A - B in directed is represented as A->B
   *        in undirected, becomes A->B and B->A
   * If the directed graph already contains reciprocal relations,
   * such as A->B and B->A, the method will throw an exception.
   */

  virtual void makeUndirected() = 0;

  // /@}


  /** @name Relations management
   *  Modificating the structure of the graph.
   */
  // /@{

  /**
   * Creates an orphaned node.
   * @return the index of the new node
   */

  virtual NodeId createNode() = 0;

  /**
   * Creates a node linked to an existing node.
   * @param origin existing node. In a directed graph: origin -> newNode.
   * @return the index of the new node
   */

  virtual NodeId createNodeFromNode(NodeId origin) = 0;

  /**
   * Creates new node on an existing Edge. A -> B will be A -> N -> B
   * @param edge existing edge.
   * @return the index of the new node
   */
  virtual NodeId createNodeOnEdge(NodeId edge) = 0;


  /**
   * Creates a node linked to new node, splitting an edge.
   * @param origin existing edge. In a directed graph: origin -> newNode.
   * @return the index of the new node
   */

  virtual NodeId createNodeFromEdge(NodeId origin) = 0;

protected:
  /**
   * Creates a link between two existing nodes. If directed graph: nodeA -> nodeB.
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   * @return the index of the new edge
   */

  virtual EdgeId link(NodeId nodeA, NodeId nodeB) = 0;

  /**
   * Sets a link between two existing nodes, using existing edge. If
   * directed graph: nodeA -> nodeB.
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   * @param edgeID the used edge
   */

  virtual void link(Graph::NodeId nodeA, Graph::NodeId nodeB, Graph::EdgeId edgeID) = 0;

  /**
   * Remove all links between two existing nodes. If directed graph: nodeA -> nodeB.
   * @param nodeA source node (or first node if undirected)
   * @param nodeB target node (or second node if undirected)
   * @return vector of IDs to de-assigned edges
   */

  virtual std::vector<EdgeId> unlink(NodeId nodeA, NodeId nodeB) = 0;

public:
  /**
   * Delete one node
   * @param node node to be deleted
   */

  virtual void deleteNode(NodeId node) = 0;

  // /@}

  /** @name Observers Management
   *  Managing communication with the observers: subscribe, unsubscribe.
   */
  // /@{

  /**
   * Attach a new observer to this Graph.
   * As a subscriber, the observer will be warned of all the changes.
   */

  virtual void registerObserver(GraphObserver* observer) = 0;

  /**
   * Detach an observer from this Graph.
   * The observer will not be warned of changes anymore.
   */

  virtual void unregisterObserver(GraphObserver* observer) = 0;

// /@}

  /** @name Nodes Functions
   *  These methodes of the graph concern the node management.
   */
  // /@{

  /**
   * @name Iterator interface on Nodes
   *
   */

  class NodeIterator
  {
public:
    virtual ~NodeIterator() {}

    virtual void next() = 0;
    virtual bool end() const = 0;
    virtual void start() = 0;

    virtual NodeId operator*() = 0;
  };

  /**
   * @brief define categories of iterators
   *
   */

  struct ALLGRAPHITER {};
  struct OUTGOINGNEIGHBORITER {};
  struct INCOMINGNEIGHBORITER {};

  /*
   * @brief builds iterator on all Nodes of the graph
   *
   */

  virtual std::unique_ptr<NodeIterator> allNodesIterator() = 0;

  virtual std::unique_ptr<NodeIterator> allNodesIterator() const = 0;

  /*
   * @brief builds iterator on the outgoing neighbor nodes of a Node
   *
   */

  virtual std::unique_ptr<NodeIterator> outgoingNeighborNodesIterator(NodeId node) = 0;

  /*
   * @brief builds iterator on the incoming neighbor nodes of a Node
   *
   */

  virtual std::unique_ptr<NodeIterator> incomingNeighborNodesIterator(NodeId node) = 0;


  /**
   * Get the number of nodes in the graph.
   */

  virtual size_t getNumberOfNodes() const = 0;

  /**
   * Get the number of edges in the graph.
   */

  virtual size_t getNumberOfEdges() const = 0;

  /**
   * Get the degree of a node (ie the number of neighbors) in the graph.
   * @param node the node one wants to count its neighbors
   * @return the number of neighbors
   */

  virtual size_t getDegree(NodeId node) const = 0;

  /**
   * Says if  a node is a leaf (ie has at most one neighbor).
   */

  virtual bool isLeaf(NodeId node) const = 0;

  /**
   * Get the number of neighbors  of a node in the graph.
   * @param node the node one wants to count its sons
   * @return the number of neighbors
   */

  virtual size_t getNumberOfNeighbors(NodeId node) const = 0;

  /**
   * Get the number of outgoing neighbors  of a node (ie the number of sons) in the graph.
   * @param node the node one wants to count its sons
   * @return the number of outgoing neighbors
   */

  virtual size_t getNumberOfOutgoingNeighbors(NodeId node) const = 0;

  /**
   * Get the number of incoming neighbors  of a node (ie the number of fathers) in the graph.
   * @param node the node one wants to count its fathers
   * @return the number of incoming neighbors
   */

  virtual size_t getNumberOfIncomingNeighbors(const NodeId node) const = 0;

  /**
   * Get all the neighbors of a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the ID of the  neighbors
   */

  virtual std::vector<NodeId> getNeighbors(const NodeId node) const = 0;

  /**
   * In an directed graph, get all the neighbors which
   * are leaving a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the ID of the  outgoing neighbors
   */

  virtual std::vector<NodeId> getOutgoingNeighbors(const NodeId node) const = 0;

  /**
   * In an directed graph, get all the neighbors which
   * are coming to a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the incoming neighbors
   */

  virtual std::vector<NodeId> getIncomingNeighbors(NodeId node) const = 0;

  /**
   * Get the leaves of a graph, ie, nodes with only one neighbor,
   * starting from a peculiar node.
   * @param node the starting node
   * @param maxDepth the maximum number of allowed depth.
   * @return a vector containing the leaves
   */

  virtual std::vector<NodeId> getLeavesFromNode(NodeId node, unsigned int maxDepth) const = 0;

  /**
   * Get all leaves of a graph, ie, nodes with no son (or only one
   * neighbor is not directet).
   * @return a vector containing the leaves
   */

  virtual std::vector<NodeId> getAllLeaves() const = 0;

  virtual std::set<NodeId> getSetOfAllLeaves() const = 0;

  /**
   * Get all the inner nodes, ie, nodes with degree > 1.
   * @return a vector containing the nodes
   */

  virtual std::vector<NodeId> getAllInnerNodes() const = 0;

  /**
   * Get all the nodes.
   * @return a vector containing the nodes
   */

  virtual std::vector<NodeId> getAllNodes() const = 0;

  /**
   * Get nodes located at the extremities of an edge
   *
   * @return a pair of the Nodes at each extremity of the edge
   *        example : N1--E1-->N2; getNodes(E1) will return (N1,N2);
   */

  virtual std::pair<NodeId, NodeId> getNodes(EdgeId edge) const = 0;

  /**
   * Get node located at the top of an edge
   *
   * @return  the Node at the top the edge
   *        example : N1--E1-->N2; getTop(E1) will return N1;
   */

  virtual NodeId getTop(EdgeId edge) const = 0;

  /**
   * Get node located at the bottom of an edge
   *
   * @return  the Node at the bottom the edge
   *        example : N1--E1-->N2; getBottom(E1) will return N2;
   */

  virtual NodeId getBottom(EdgeId edge) const = 0;

  // /@}

  /** @name Topological Properties
   *  These methodes check some topological properties.
   */
  // /@{

  /**
   * Is the graph a tree?
   * @return false if a node is met more than one time browsing the graph
   */

  virtual bool isTree() const = 0;

  /**
   * Is the graph directed acyclic?
   * @return true if an edge is met more than one time browsing the graph
   */

  virtual bool isDA() const = 0;


  /**
   * Orientates the graph hanging from the root
   */

  virtual void orientate() = 0;

  /**
   * Is the graph directed?
   * @return true the type of the graph is directed
   */

  virtual bool isDirected() const = 0;


  /**
   * Does the graph contain reciprocal relations such as A->B and B->A?
   * @return true if one of them is seen in the structure
   */

  virtual bool containsReciprocalRelations() const = 0;


  // /@}


  /** @name Edge Functions
   *  These methodes of the graph concern the edges.
   */
  // /@{

  /**
   * @name Iterator interface on Nodes
   *
   */

  class EdgeIterator
  {
public:
    virtual ~EdgeIterator() {}

    virtual void next() = 0;
    virtual bool end() const = 0;
    virtual void start() = 0;

    virtual EdgeId operator*() = 0;
  };


  /*
   * @brief builds iterator on all Nodes of the graph
   *
   */

  virtual std::unique_ptr<EdgeIterator> allEdgesIterator() = 0;

  /*
   * @brief builds iterator on the outgoing neighbor nodes of a Node
   *
   */

  virtual std::unique_ptr<EdgeIterator> outgoingEdgesIterator(NodeId node) = 0;

  /*
   * @brief builds iterator on the incoming neighbor nodes of a Node
   *
   */

  virtual std::unique_ptr<EdgeIterator> incomingEdgesIterator(NodeId node) = 0;

  /**
   * Get all the edges to/from a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the ID of the  edges
   */

  virtual std::vector<EdgeId> getEdges(const NodeId node) const = 0;

  /**
   * In an directed graph, get all the edges which
   * are leaving a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the ID of the outgoing edges
   */

  virtual std::vector<EdgeId> getOutgoingEdges(const NodeId node) const = 0;

  /**
   * In an directed graph, get all the edges which
   * are coming to a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the incoming edges
   */

  virtual std::vector<EdgeId> getIncomingEdges(NodeId node) const = 0;

  /**
   * Get all edges of a graph.
   * @return a vector containing the edges
   */

  virtual std::vector<EdgeId> getAllEdges() const = 0;

  /**
   * Returns the Edge between two nodes
   * @param nodeA if directed, origin node
   * @param nodeB if directed, destination node
   * @return the edge between these two nodes
   */

  virtual EdgeId getEdge(NodeId nodeA, NodeId nodeB) const = 0;

  /**
   * Returns the Edge between two nodes, trying both directions
   * @param nodeA any node implied in the relation
   * @param nodeB any other node implied in the relation
   * @return the edge between these two nodes
   */

  virtual EdgeId getAnyEdge(NodeId nodeA, NodeId nodeB) const = 0;

  // /@}

protected:
  /** @name Updating the changes on the observers
   *  These methodes aim to trigger some changes to the observers
   */

  // /@{

  /**
   * Trigger E objects deleting on the observers
   * @param edgesToDelete list of edges to delete
   */
  virtual void notifyDeletedEdges(const std::vector<EdgeId>& edgesToDelete) const = 0;

  /**
   * Trigger N objects deleting on the observers
   * @param nodesToDelete list of edges to delete
   */
  virtual void notifyDeletedNodes(const std::vector<NodeId>& nodesToDelete) const = 0;


  // /@}

public:
  /**
   * Output the graph in DOT format
   * @param out a ostream where the DOT format will be output
   * @param name a string naming the graph
   */

  virtual void outputToDot(std::ostream& out, const std::string& name) const = 0;

  template<class N, class E, class GraphImpl>
  friend class AssociationGraphImplObserver;
};
}

#endif
