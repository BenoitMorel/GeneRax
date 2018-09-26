//
// File AssociationGraphObserver.h
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

#ifndef _ASSOCIATIONGRAPHOBSERVER_HPP_
#define _ASSOCIATIONGRAPHOBSERVER_HPP_

#include "GraphObserver.h"

#include "../Exceptions.h"
#include "../Clonable.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <memory>
#include <type_traits>

namespace bpp
{
/**
 * @brief Defines a Graph Associator. It is a template which follows
 * (subscribed to) a Graph.
 *
 * @author Thomas Bigot, Laurent Guéguen
 */

template<class N, class E>
class AssociationGraphObserver :
  public virtual GraphObserver
{
private:
  typedef Graph::NodeId NodeGraphid;
  typedef Graph::EdgeId EdgeGraphid;

public:
  typedef unsigned int NodeIndex;
  typedef unsigned int EdgeIndex;

  /**
   * Specific copy of A and B objects, if clonable or not.
   *
   */

  template<class A, class B>
  static B* copy(const A& a, typename std::enable_if< !std::is_base_of<B, A>::value&& !std::is_abstract<B>::value>::type* = 0)
  {
    return new B(a);
  }

  template<class A, class B>
  static B* copy(const A& a, typename std::enable_if< std::is_base_of<B, A>::value&& !std::is_abstract<A>::value>::type* = 0)
  {
    return dynamic_cast<B*>(new A(a));
  }

  template<class A, class B>
  static B* copy(const A& a, typename std::enable_if< std::is_base_of<B, A>::value&& std::is_abstract<A>::value&& std::is_base_of<Clonable, A>::value>::type* = 0)
  {
    return dynamic_cast<B*>(a.clone());
  }


  /** @name Graph Relations Management
   *  Modificating the structure of the graph.
   */
  // /@{

  /**
   * Creates an orphaned node from a NodeClass object.
   * @param newNodeObject the N object associated to the node in the graph.
   *
   */
  virtual void createNode(std::shared_ptr<N>  newNodeObject) = 0;


  /**
   * Creates an node linked to an existing node. Order of parameters match
   * the link method.
   * @param newNodeObject the N object associated to the node in the graph.
   * @param objectOriginNode existing node. In a directed graph:
   * origin -> newNode.
   * @param newEdgeObject optional edge between nodes (default = 00)
   */

  virtual void createNode(std::shared_ptr<N>  objectOriginNode, std::shared_ptr<N>  newNodeObject, std::shared_ptr<E>  newEdgeObject = 00) = 0;

  /**
   * Creates a link between two existing nodes.
   * If directed graph: nodeA -> nodeB.
   * @param nodeObjectA source node (or first node if undirected)
   * @param nodeObjectB target node (or second node if undirected)
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */

  virtual void link(std::shared_ptr<N>  nodeObjectA, std::shared_ptr<N>  nodeObjectB, std::shared_ptr<E>  edgeObject = 00) = 0;

  /**
   * Detroys a link between two existing nodes in the graph.
   * If directed graph: nodeA -> nodeB.
   * @param nodeObjectA source node (or first node if undirected)
   * @param nodeObjectB target node (or second node if undirected)
   */

  virtual void unlink(std::shared_ptr<N>  nodeObjectA, std::shared_ptr<N>  nodeObjectB) = 0;

  /**
   * Deletes a node
   * @param nodeObject node to be deleted
   */
  virtual void deleteNode(std::shared_ptr<N>  nodeObject) = 0;

  // /@}

  /** @name Object Association
   *  Associate or dissociate N and E objects to pre-existing Graph Nodes and Graph Edges
   */
  // /@{

  /**
   * Associate a N object to a node in the graph.
   * @param nodeObject object to associate
   * @param node existing node to be associated
   */

  virtual void associateNode(std::shared_ptr<N>  nodeObject, NodeGraphid node) = 0;

  /**
   * Associate a E object to an edge in the graph.
   * @param edgeObject object to associate
   * @param edge existing edge to be associated
   */

  virtual void associateEdge(std::shared_ptr<E>  edgeObject, EdgeGraphid edge) = 0;

  /**
   * Dissociate a N or E object to a node or an edge in the graph.
   * Graph is not changed.
   *
   * @param nodeObject object to dissociate
   */

  virtual void dissociateNode(std::shared_ptr<N>  nodeObject) = 0;
  virtual void dissociateEdge(std::shared_ptr<E>  edgeObject) = 0;


  /**
   * Return the associated Node ID
   * @param nodeObject object which to return the node ID
   * @return a node ID
   */
  virtual NodeGraphid getNodeGraphid(const std::shared_ptr<N>  nodeObject) const = 0;

  /**
   * Return the associated Node ID
   * @param edgeObject object which to return the node ID
   * @return a node ID
   */
  virtual EdgeGraphid getEdgeGraphid(const std::shared_ptr<E>  edgeObject) const = 0;

  /**
   * Transforms an (a list of) id(s) into an (a list of) object(s)
   */
  virtual const std::shared_ptr<N>  getNodeFromGraphid(NodeGraphid) const = 0;
  virtual std::shared_ptr<N>  getNodeFromGraphid(NodeGraphid) = 0;

  virtual std::vector<std::shared_ptr<N> > getNodesFromGraphid(std::vector<NodeGraphid> ) const = 0;
  virtual std::shared_ptr<E>  getEdgeFromGraphid(EdgeGraphid) = 0;
  virtual const std::shared_ptr<E>  getEdgeFromGraphid(EdgeGraphid) const = 0;
  virtual std::vector<std::shared_ptr<E> > getEdgesFromGraphid(std::vector<EdgeGraphid> ) const = 0;


  /**
   * @return the root
   */

  virtual std::shared_ptr<N> getRoot() const = 0;

  virtual NodeIndex getRootIndex() const = 0;

  // /@}


  /** @name Object Indexation
   *  Get or set indexes to nodes and edges
   */
  // /@{

  /**
   * @brief return if the object has an index.
   */

  virtual bool hasIndex(const std::shared_ptr<N> nodeObject) const = 0;
  virtual bool hasIndex(const std::shared_ptr<E> edgeObject) const = 0;

  /**
   * Return the associated Node index
   * @param nodeObject object which to return the node index
   * @return a node index
   */

  virtual NodeIndex getNodeIndex(const std::shared_ptr<N>  nodeObject) const = 0;
  virtual std::vector<NodeIndex> getNodeIndexes(std::vector<std::shared_ptr<N> > nodeObjects) const = 0;


  /**
   * Return the associated Node index
   * @param edgeObject object which to return the node index
   * @return a node index
   */
  virtual EdgeIndex getEdgeIndex(const std::shared_ptr<E>  edgeObject) const = 0;
  virtual std::vector<EdgeIndex> getEdgeIndexes(std::vector<std::shared_ptr<E> > edgeObjects) const = 0;

  /**
   * Set an index associated to a node
   * @param nodeObject object to which one want to set the index
   * @param index intex to be given, 0 to get the first free index
   * @return the given index
   */
  virtual NodeIndex setNodeIndex(const std::shared_ptr<N>  nodeObject, NodeIndex index) = 0;

  /**
   * Set an index associated to an edge
   * @param edgeObject object to which one want to set the index
   * @param index intex to be given, 0 to get the first free index
   * @return the given index
   */
  virtual EdgeIndex setEdgeIndex(const std::shared_ptr<E>  edgeObject, EdgeIndex index) = 0;

  /**
   * Return the associated Node, querying with an index
   * @param nodeIndex the index of the wanted node
   * @return N, a node object
   */

  virtual std::shared_ptr<N>  getNode(NodeIndex nodeIndex) const = 0;

  /**
   * Return the associated Node index
   * @param edgeIndex the index of the wanted edge
   * @return E, an edge object
   */
  virtual std::shared_ptr<E>  getEdge(EdgeIndex edgeIndex) const = 0;

  // /@}

  /** @name Topology exploration
   *  These methodes of the graph concern the topology exploration.
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

    virtual std::shared_ptr<N> operator*() = 0;
  };


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

  virtual std::unique_ptr<NodeIterator> outgoingNeighborNodesIterator(std::shared_ptr<N> node) = 0;

  virtual std::unique_ptr<NodeIterator> outgoingNeighborNodesIterator(std::shared_ptr<N> node) const = 0;

  /*
   * @brief builds iterator on the incoming neighbor nodes of a Node
   *
   */

  virtual std::unique_ptr<NodeIterator> incomingNeighborNodesIterator(std::shared_ptr<N> node) = 0;

  virtual std::unique_ptr<NodeIterator> incomingNeighborNodesIterator(std::shared_ptr<N> node) const = 0;

  /**
   * Get all the neighbors of a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the neighbors
   */

  virtual std::vector<std::shared_ptr<N> > getNeighbors(const std::shared_ptr<N>  node) const = 0;
  virtual std::vector<NodeIndex> getNeighbors(NodeIndex node) const = 0;

  /**
   * Get all the edges from/to a node in the graph.
   * @param node the node one wants to get its neighbor edges
   * @return a vector containing the edges
   */
  virtual std::vector<std::shared_ptr<E> > getEdges(const std::shared_ptr<N>  node) const = 0;
  virtual std::vector<EdgeIndex> getEdges(NodeIndex node) const = 0;

  /**
   * In an directed graph, get all the neighbors which
   * are leaving a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the outgoing neighbors
   */
  virtual std::vector<std::shared_ptr<N> > getOutgoingNeighbors(const std::shared_ptr<N>  node) const = 0;
  virtual std::vector<NodeIndex> getOutgoingNeighbors(NodeIndex node) const = 0;

  /**
   * In an directed graph, get all the edges which
   * are leaving a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the outgoing edges
   */
  virtual std::vector<std::shared_ptr<E> > getOutgoingEdges(const std::shared_ptr<N>  node) const = 0;
  virtual std::vector<EdgeIndex> getOutgoingEdges(NodeIndex node) const = 0;


  /**
   * In an directed graph, get all the neighbors which
   * are coming to a node in the graph.
   * @param node the node one wants to get its neighbors
   * @return a vector containing the incoming neighbors
   */
  virtual std::vector<std::shared_ptr<N> > getIncomingNeighbors(const std::shared_ptr<N>  node) const = 0;
  virtual std::vector<NodeIndex> getIncomingNeighbors(NodeIndex node) const = 0;

  /**
   * In an directed graph, get all the edges which
   * are coming to a node in the graph.
   * @param node the node one wants to get its edges
   * @return a vector containing the incoming edges
   */
  virtual std::vector<std::shared_ptr<E> > getIncomingEdges(const std::shared_ptr<N>  node) const = 0;
  virtual std::vector<EdgeIndex> getIncomingEdges(NodeIndex node) const = 0;


  /**
   * Get the leaves of a graph, ie, nodes with only one neighbor,
   * starting from a peculiar node, up to a specific depth.
   * @param node the starting node
   * @param maxDepth the max recursion depth
   * @return a vector containing the leaves
   */
  virtual std::vector<std::shared_ptr<N> > getLeavesFromNode(std::shared_ptr<N>  node, unsigned int maxDepth) const = 0;

  /**
   * Get all the leaves objects of a graph, ie, nodes with only one neighbor,
   * @return a vector containing the leaves
   */
  virtual std::vector<std::shared_ptr<N> > getAllLeaves() const = 0;

  /**
   * Get all the inner nodes of a graph, ie, nodes with degree >= 2.
   * @return a vector containing the inner nodes.
   */
  virtual std::vector<std::shared_ptr<N> > getAllInnerNodes() const = 0;

  /**
   * Get all the defined nodes of a graphO,
   * @return a vector containing the nodesObjects
   */
  virtual std::vector<std::shared_ptr<N> > getAllNodes() const = 0;
  virtual std::vector<NodeIndex> getAllNodesIndexes() const = 0;

  /**
   * @name Iterator interface on Edges
   *
   */

  class EdgeIterator
  {
public:
    virtual ~EdgeIterator() {}

    virtual void next() = 0;
    virtual bool end() const = 0;
    virtual void start() = 0;

    virtual std::shared_ptr<E> operator*() = 0;
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

  virtual std::unique_ptr<EdgeIterator> outgoingEdgesIterator(std::shared_ptr<N> node) = 0;

  /*
   * @brief builds iterator on the incoming neighbor nodes of a Node
   *
   */

  virtual std::unique_ptr<EdgeIterator> incomingEdgesIterator(std::shared_ptr<N> node) = 0;

  /**
   * Get all the defined edges of a graphO,
   * @return a vector containing the branchesObjects
   */

  virtual std::vector<std::shared_ptr<E> > getAllEdges() const = 0;

  /**
   * Returns the Edge between two nodes nodeA -> nodeB
   * @param nodeA source node (if directed)
   * @param nodeB destination node (if directed)
   * @return the edge between these two nodes
   */
  virtual std::shared_ptr<E>  getEdgeLinking(std::shared_ptr<N>  nodeA, std::shared_ptr<N>  nodeB) const = 0;

  /**
   * Sets the Edge between two nodes nodeA -> nodeB
   * @param nodeA source node (if directed)
   * @param nodeB destination node (if directed)
   * @param edge the edge between these two nodes
   */
  virtual void setEdgeLinking(std::shared_ptr<N>  nodeA, std::shared_ptr<N>  nodeB, std::shared_ptr<E>  edge) = 0;

  // /@}
};
}

#endif
