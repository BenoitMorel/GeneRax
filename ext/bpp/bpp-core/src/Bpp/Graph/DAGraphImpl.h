//
// File DAGraphImpl.h
// Created by: Laurent Guéguen
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

#ifndef _DA_GRAPH_IMPL_H_
#define _DA_GRAPH_IMPL_H_

#include <string>
#include <vector>
#include <iostream>
#include <ostream>


#include "DAGraph.h"
#include "GlobalGraph.h"

#include "../Exceptions.h"
#include "../Numeric/VectorTools.h"

namespace bpp
{
template<class GraphImpl>
class DAGraphImpl :
  public virtual DAGraph,
  public GraphImpl
{
protected:
  /**
   * Is the graph a DAG? Set to false when structure is modified, true after validation.
   */

  mutable bool isValid_;

  /**
   * Is the graph rooted? Set to false when structure is modified,
   * true after validation.
   */

  mutable bool isRooted_;

  // unvalidate the DAG
  virtual void topologyHasChanged_() const;

  // will throw an exception if the DA is not valid
  void mustBeValid_() const;

  // will throw an exception if the DA is not rooted, ie as more
  // than one node with no father.

  void mustBeRooted_() const;

  // test the validity of the DAG
  bool validate_() const;

  /**
   * Reorient at mimina all the edges starting from a node: the
   * father nodes become sons, and so on.
   *
   * The Graph must already be rooted
   */

  void propagateDirection_(Graph::NodeId node);

  /**
   * Reorient all the graph so that localRoot is ahead the other
   * sons, with the exception of the nodes already met
   *
   */

  void orientGraphFrom_(std::set<Graph::NodeId>& metNodes, Graph::NodeId localRoot);

public:
  /**
   * bool is only for inheritance from observers, useless.
   *
   */

  DAGraphImpl(bool b = true);

  /**
   * Is the graph a DAG?
   * @return true if valid DAG
   */
  bool isValid() const;

  /**
   * Is the DAG rooted?
   *
   * @return true if rooted, ie has only one node with no father.
   */

  bool isRooted() const;

  /**
   * Says if  a node is a leaf (ie has at most one neighbor).
   */

  bool isLeaf(Graph::NodeId node) const;

  /**
   * Check if node has a father
   */

  bool hasFather(Graph::NodeId node) const;

  /**
   * Get the father nodes of a node
   * @return the father node
   */

  std::vector<Graph::NodeId> getFathers(Graph::NodeId nodeid) const;

  /**
   * @brief Get the number of fathers nodes
   */

  size_t getNumberOfFathers(Graph::NodeId node) const;

  /**
   * Add a father to a node
   */

  void addFather(Graph::NodeId node, Graph::NodeId father);

  void addFather(Graph::NodeId node, Graph::NodeId father, Graph::EdgeId edgeId);

  /**
   * Remove one father
   */

  void removeFather(Graph::NodeId node, Graph::NodeId father);

  /**
   * Remove all the fathers
   */

  std::vector<Graph::NodeId> removeFathers(Graph::NodeId node);

  /**
   * Get the leaves under a node
   * @param node the starting node
   * @return a vector containing the leaves
   */

  std::vector<Graph::NodeId> getLeavesUnderNode(Graph::NodeId node) const;

  /**
   * Get the sons node of a node
   */

  std::vector<Graph::NodeId> getSons(Graph::NodeId node) const;

  /**
   * @brief Get the number of sons node
   */

  size_t getNumberOfSons(Graph::NodeId node) const;

  /**
   * Add a son to a node
   */

  void addSon(Graph::NodeId node, Graph::NodeId sonNode);

  void addSon(Graph::NodeId node, Graph::NodeId sonNode, Graph::EdgeId edge);

  /**
   * Remove all the sons
   */

  std::vector<Graph::NodeId> removeSons(Graph::NodeId node);

  /**
   * Remove one son
   */

  void removeSon(Graph::NodeId node, Graph::NodeId son);

  /**
   * Re-root the DA with the new root (and make the graph a DA if it is not)
   */

  void rootAt(Graph::NodeId newRoot);

  /**
   * Get all the nodes below a node
   */

  std::vector<Graph::NodeId> getBelowNodes(Graph::NodeId localRoot) const;

  /**
   * Get all the branches below a node
   */

  std::vector<Graph::EdgeId> getBelowEdges(Graph::NodeId localRoot) const;

  // recursive function for getSubtreeNodes
  void fillSubtreeMetNodes_(std::vector<Graph::NodeId>& metNodes, Graph::NodeId localRoot) const;

  // recursive function for getSubtreeEdges
  void fillSubtreeMetEdges_(std::vector<Graph::EdgeId>& metEdges, Graph::NodeId localRoot) const;

  // recursive function for getLeavesUnderNode
  void fillListOfLeaves_(Graph::NodeId startingNode, std::vector<Graph::NodeId>& foundLeaves) const;
};


/******************/

typedef DAGraphImpl<GlobalGraph> DAGlobalGraph;

/*****************/


template<class GraphImpl>
DAGraphImpl<GraphImpl>::DAGraphImpl(bool b) :
  GraphImpl(true),
  isValid_(false),
  isRooted_(false)
{}


template<class GraphImpl>
bool DAGraphImpl<GraphImpl>::isValid() const
{
  return isValid_ || validate_();
}


template<class GraphImpl>
bool DAGraphImpl<GraphImpl>::isRooted() const
{
  if (isRooted_)
    return true;

  std::unique_ptr<Graph::NodeIterator> allIt = allNodesIterator();
  bool seen = false;

  for ( ; !allIt->end(); allIt->next())
  {
    if (getNumberOfFathers(**allIt) == 0)
    {
      if (seen)
        return false;
      seen = true;
    }
  }

  isRooted_ = seen;
  return true;
}

template<class GraphImpl>
bool DAGraphImpl<GraphImpl>::isLeaf(Graph::NodeId node) const
{
  return GraphImpl::isLeaf(node) == 0;
}

template<class GraphImpl>
bool DAGraphImpl<GraphImpl>::hasFather(Graph::NodeId node) const
{
  return GraphImpl::getNumberOfIncomingNeighbors(node) >= 1;
}

template<class GraphImpl>
std::vector<Graph::NodeId> DAGraphImpl<GraphImpl>::getFathers(Graph::NodeId node) const
{
  return getIncomingNeighbors(node);
}

template<class GraphImpl>
size_t DAGraphImpl<GraphImpl>::getNumberOfFathers(Graph::NodeId node) const
{
  return GraphImpl::getNumberOfIncomingNeighbors(node);
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::addFather(Graph::NodeId node, Graph::NodeId fatherNode)
{
  GraphImpl::link(fatherNode, node);
  topologyHasChanged_();
  isRooted_ = false;
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::addFather(Graph::NodeId node, Graph::NodeId fatherNode, Graph::EdgeId edge)
{
  GraphImpl::link(fatherNode, node, edge);
  topologyHasChanged_();
  isRooted_ = false;
}

template<class GraphImpl>
std::vector<Graph::NodeId> DAGraphImpl<GraphImpl>::removeFathers(Graph::NodeId node)
{
  std::vector<Graph::NodeId> fathers = getFathers(node);
  for (std::vector<Graph::NodeId>::iterator currFather = fathers.begin(); currFather != fathers.end(); currFather++)
  {
    removeFather(node, *currFather);
  }
  return fathers;
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::removeFather(Graph::NodeId node, Graph::NodeId father)
{
  if (getNumberOfIncomingNeighbors(node) == 1)
    isRooted_ = false;

  GraphImpl::unlink(father, node);
}

template<class GraphImpl>
std::vector<Graph::NodeId> DAGraphImpl<GraphImpl>::getLeavesUnderNode(Graph::NodeId node) const
{
  std::vector<Graph::NodeId> foundLeaves;
  fillListOfLeaves_(node, foundLeaves);

  return foundLeaves;
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::fillListOfLeaves_(Graph::NodeId startingNode, std::vector<Graph::NodeId>& foundLeaves) const
{
  const std::vector<Graph::NodeId> sons = getSons(startingNode);
  if (sons.size() > 1)
  {
    for (std::vector<Graph::NodeId>::const_iterator currNeighbor = sons.begin(); currNeighbor != sons.end(); currNeighbor++)
    {
      fillListOfLeaves_(*currNeighbor, foundLeaves);
    }
  }
  else
  {
    foundLeaves.push_back(startingNode);
  }
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::mustBeRooted_() const
{
  if (!isRooted())
    throw Exception("DAGraphImpl<GraphImpl>: The DAG must be rooted.");
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::mustBeValid_() const
{
  if (!isValid())
    throw Exception("DAGraphImpl<GraphImpl>: The DAG is not valid.");
}


template<class GraphImpl>
bool DAGraphImpl<GraphImpl>::validate_() const
{
  isValid_ = GraphImpl::isDA();
  return isValid_;
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::topologyHasChanged_() const
{
  isValid_ = false;
}


template<class GraphImpl>
std::vector<Graph::NodeId> DAGraphImpl<GraphImpl>::getSons(Graph::NodeId node) const
{
  return GraphImpl::getOutgoingNeighbors(node);
}

template<class GraphImpl>
size_t DAGraphImpl<GraphImpl>::getNumberOfSons(Graph::NodeId node) const
{
  return GraphImpl::getNumberOfOutgoingNeighbors(node);
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::addSon(Graph::NodeId node, Graph::NodeId sonNode)
{
  GraphImpl::link(node, sonNode);
  topologyHasChanged_();
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::addSon(Graph::NodeId node, Graph::NodeId sonNode, Graph::EdgeId edge)
{
  GraphImpl::link(node, sonNode, edge);
  topologyHasChanged_();
}


template<class GraphImpl>
std::vector<Graph::NodeId> DAGraphImpl<GraphImpl>::removeSons(Graph::NodeId node)
{
  std::vector<Graph::NodeId> sons = getSons(node);
  for (std::vector<Graph::NodeId>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
  {
    removeSon(node, *currSon);
  }
  return sons;
}


template<class GraphImpl>
void DAGraphImpl<GraphImpl>::removeSon(Graph::NodeId node, Graph::NodeId son)
{
  GraphImpl::unlink(node, son);
}


template<class GraphImpl>
void DAGraphImpl<GraphImpl>::rootAt(Graph::NodeId newRoot)
{
  GraphImpl::setRoot(newRoot);

  // change edge direction between the new node and the former fathers
  if (isRooted() && isValid())
    propagateDirection_(newRoot);
  else
  {
    GraphImpl::orientate();
    isRooted_ = true;
  }
}


template<class GraphImpl>
void DAGraphImpl<GraphImpl>::propagateDirection_(Graph::NodeId node)
{
  std::vector<Graph::NodeId> vFat = getFathers(node);

  for (size_t i = 0; i < vFat.size(); i++)
  {
    propagateDirection_(vFat[i]);
  }

  for (size_t i = 0; i < vFat.size(); i++)
  {
    GraphImpl::switchNodes(vFat[i], node);
  }
}


template<class GraphImpl>
std::vector<Graph::NodeId> DAGraphImpl<GraphImpl>::getBelowNodes(Graph::NodeId localRoot) const
{
  mustBeValid_();
  std::vector<Graph::EdgeId> metNodes;
  fillSubtreeMetNodes_(metNodes, localRoot);
  return metNodes;
}

template<class GraphImpl>
std::vector<Graph::EdgeId> DAGraphImpl<GraphImpl>::getBelowEdges(Graph::NodeId localRoot) const
{
  mustBeValid_();
  std::vector<Graph::EdgeId> metEdges;
  fillSubtreeMetEdges_(metEdges, localRoot);
  return metEdges;
}


template<class GraphImpl>
void DAGraphImpl<GraphImpl>::fillSubtreeMetNodes_(std::vector<Graph::NodeId>& metNodes, Graph::NodeId localRoot) const
{
  metNodes.push_back(localRoot);
  std::vector<Graph::NodeId> sons = GraphImpl::getOutgoingNeighbors(localRoot);
  for (std::vector<Graph::NodeId>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
  {
    fillSubtreeMetNodes_(metNodes, *currSon);
  }
}

template<class GraphImpl>
void DAGraphImpl<GraphImpl>::fillSubtreeMetEdges_(std::vector<Graph::EdgeId>& metEdges, Graph::NodeId localRoot) const
{
  std::vector<Graph::EdgeId> edgesToSons = GraphImpl::getOutgoingEdges(localRoot);
  for (std::vector<Graph::EdgeId>::iterator currEdgeToSon = edgesToSons.begin(); currEdgeToSon != edgesToSons.end(); currEdgeToSon++)
  {
    metEdges.push_back(*currEdgeToSon);
    fillSubtreeMetEdges_(metEdges, GraphImpl::getBottom(*currEdgeToSon));
  }
}
}


#endif
