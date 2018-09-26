//
// File TreeGraphImpl.h
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

#ifndef _TREEGRAPH_IMPL_H_
#define _TREEGRAPH_IMPL_H_

#include <string>
#include <vector>
#include <iostream>
#include <ostream>


#include "TreeGraph.h"
#include "GlobalGraph.h"

#include "../Exceptions.h"
#include "../Numeric/VectorTools.h"

namespace bpp
{
template<class GraphImpl>
class TreeGraphImpl :
  public virtual TreeGraph,
  public GraphImpl
{
private:
  /**
   * Is the graph a tree? Set to false when structure is modified, true after validation.
   */
  mutable bool isValid_;

  // unvalidate the tree
  void topologyHasChanged_() const;

  // will throw an exception if the tree is not valid
  void mustBeValid_() const;

  // will throw an exception if the tree is not rooted
  void mustBeRooted_() const;

  // test the validity of the tree
  bool validate_() const;

  /**
   * Reorient all the edges starting from a node:
   * the father node becomes a son, and so on.
   */
  void propagateDirection_(Graph::NodeId node);

  // recursive function for getSubtreeNodes
  void fillSubtreeMetNodes_(std::vector<Graph::NodeId>& metNodes, Graph::NodeId localRoot) const;

  // recursive function for getSubtreeEdges
  void fillSubtreeMetEdges_(std::vector<Graph::EdgeId>& metEdges, Graph::NodeId localRoot) const;

  // recursive function for getLeavesUnderNode
  void fillListOfLeaves_(Graph::NodeId startingNode, std::vector<Graph::NodeId>& foundLeaves) const;

public:
  TreeGraphImpl();

  TreeGraphImpl(bool rooted = true);

  /**
   * Is the graph a tree? A tree must be acyclic and with no isolated node.
   * @return true if valid tree
   */
  bool isValid() const;

  /**
   * Is the tree rooted?
   *
   * @return true if rooted, ie is directed.
   */

  bool isRooted() const;

  /**
   * Get the father node of a node in a rooted tree
   * @return the father node
   */

  Graph::NodeId getFather(Graph::NodeId nodeid) const;

  /**
   * Get the branch leading to the father in a rooted tree
   * @return the branch between a node and its father
   */

  Graph::EdgeId getEdgeToFather(Graph::NodeId node) const;

  /**
   * Check if node has a father
   */

  bool hasFather(Graph::NodeId node) const;

  /**
   * Says if  a node is a leaf (ie has at most one neighbor).
   */

  bool isLeaf(Graph::NodeId node) const;

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
   * Get the branches to the sons node of a node
   */

  std::vector<Graph::EdgeId> getBranches(Graph::NodeId node) const;

  /**
   * Get a iterator on the sons node of a node
   */

  std::unique_ptr<Graph::NodeIterator> sonsIterator(Graph::NodeId node);

  /**
   * Get a iterator on the branches to sons of a node
   */

  std::unique_ptr<Graph::EdgeIterator> branchesIterator(Graph::NodeId node);

  /**
   * @brief Get the number of sons node
   */

  size_t getNumberOfSons(Graph::NodeId node) const;

  /**
   * set the father node of a node in a rooted tree
   */

  void setFather(Graph::NodeId node, Graph::NodeId fatherNode);
  void setFather(Graph::NodeId node, Graph::NodeId fatherNode, Graph::EdgeId edgeId);

  /**
   * Add a son to a node in a rooted tree
   */

  void addSon(Graph::NodeId node, Graph::NodeId sonNode);

  void addSon(Graph::NodeId node, Graph::NodeId sonNode, Graph::EdgeId edgeId);

  /**
   * Remove all the sons
   */

  std::vector<Graph::NodeId> removeSons(Graph::NodeId node);

  /**
   * Remove one son
   */

  void removeSon(Graph::NodeId node, Graph::NodeId son);

  /**
   * Re-root the tree with the new root
   */

  void rootAt(Graph::NodeId newRoot);

  /**
   * Set the tree to its flat unrooted version.
   * As an algorithmical convenience, a root node is kept, but it has
   * no logical significance.
   */

  void unRoot(bool joinRootSons);

  /**
   * Set a node as a new outgroup in a rooted tree, will make a root between
   * the given node and its father.
   */

  void setOutGroup(Graph::NodeId newOutGroup);

  /**
   * Get all the nodes of a subtree
   */

  std::vector<Graph::NodeId> getSubtreeNodes(Graph::NodeId localRoot) const;

  /**
   * Get all the branches of a subtree
   */

  std::vector<Graph::EdgeId> getSubtreeEdges(Graph::NodeId localRoot) const;

  // ///FROM TREETOOLS & TREETOOLS COMPAT


  std::vector<Graph::NodeId> getNodePathBetweenTwoNodes(Graph::NodeId nodeA, Graph::NodeId nodeB, bool includeAncestor = true) const;
  std::vector<Graph::EdgeId> getEdgePathBetweenTwoNodes(Graph::NodeId nodeA, Graph::NodeId nodeB) const;
};


/******************/

typedef TreeGraphImpl<GlobalGraph> TreeGlobalGraph;

/*****************/


template<class GraphImpl>
TreeGraphImpl<GraphImpl>::TreeGraphImpl(bool rooted) :
  GraphImpl(rooted),
  isValid_(false)
{}


template<class GraphImpl>
bool TreeGraphImpl<GraphImpl>::isValid() const
{
  return isValid_ || validate_();
}

template<class GraphImpl>
Graph::NodeId TreeGraphImpl<GraphImpl>::getFather(Graph::NodeId node) const
{
  std::vector<Graph::NodeId> incomers = getIncomingNeighbors(node);
  if (incomers.size() > 1)
    throw Exception("TreeGraphImpl<GraphImpl>::getFather: more than one father for Node " + TextTools::toString(node) + " : " + VectorTools::paste(incomers, ",") + ". Should never happen since validity has been controled. Please report this bug.");
  if (incomers.size() == 0)
    throw Exception("TreeGraphImpl<GraphImpl>::getFather: node " + TextTools::toString(node) + " has no father.");
  return *incomers.begin();
}

template<class GraphImpl>
Graph::EdgeId TreeGraphImpl<GraphImpl>::getEdgeToFather(Graph::NodeId node) const
{
  Graph::NodeId father = getFather(node);
  return GraphImpl::getEdge(father, node);
}

template<class GraphImpl>
bool TreeGraphImpl<GraphImpl>::hasFather(Graph::NodeId node) const
{
  return GraphImpl::getNumberOfIncomingNeighbors(node) >= 1;
}

template<class GraphImpl>
bool TreeGraphImpl<GraphImpl>::isLeaf(Graph::NodeId node) const
{
  return (!GraphImpl::isDirected() && GraphImpl::getNumberOfOutgoingNeighbors(node) <= 1)
         || (GraphImpl::isDirected() && GraphImpl::getNumberOfOutgoingNeighbors(node) == 0);
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::fillListOfLeaves_(Graph::NodeId startingNode, std::vector<Graph::NodeId>& foundLeaves) const
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
std::vector<Graph::NodeId> TreeGraphImpl<GraphImpl>::getLeavesUnderNode(Graph::NodeId node) const
{
  std::vector<Graph::NodeId> foundLeaves;
  fillListOfLeaves_(node, foundLeaves);

  return foundLeaves;
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::mustBeRooted_() const
{
  if (!isRooted())
    throw Exception("TreeGraphImpl<GraphImpl>: The tree must be rooted.");
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::mustBeValid_() const
{
  if (!isValid())
    throw Exception("TreeGraphImpl<GraphImpl>: The tree is not valid.");
}

template<class GraphImpl>
bool TreeGraphImpl<GraphImpl>::isRooted() const
{
  return GraphImpl::isDirected();
}

template<class GraphImpl>
bool TreeGraphImpl<GraphImpl>::validate_() const
{
  isValid_ = GraphImpl::isTree();
  return isValid_;
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::topologyHasChanged_() const
{
  isValid_ = false;
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::rootAt(Graph::NodeId newRoot)
{
  if (!isValid())
    throw Exception("TreeGraphImpl::rootAt: Tree is not Valid.");

  GraphImpl::makeDirected();
  // set the new root on the Graph
  GraphImpl::setRoot(newRoot);
  // change edge direction between the new node and the former one
  propagateDirection_(newRoot);
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::propagateDirection_(Graph::NodeId node)
{
  if (hasFather(node))
  {
    NodeId father = getFather(node);
    propagateDirection_(father);
    GraphImpl::switchNodes(father, node);
  }
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::setFather(Graph::NodeId node, Graph::NodeId fatherNode)
{
  if (hasFather(node))
    GraphImpl::unlink(getFather(node), node);
  GraphImpl::link(fatherNode, node);
  topologyHasChanged_();
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::setFather(Graph::NodeId node, Graph::NodeId fatherNode, Graph::EdgeId edgeId)
{
  if (hasFather(node))
    GraphImpl::unlink(getFather(node), node);
  GraphImpl::link(fatherNode, node, edgeId);
  topologyHasChanged_();
}


template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::addSon(Graph::NodeId node, Graph::NodeId sonNode)
{
  GraphImpl::link(node, sonNode);
  topologyHasChanged_();
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::addSon(Graph::NodeId node, Graph::NodeId sonNode, Graph::EdgeId edgeId)
{
  GraphImpl::link(node, sonNode, edgeId);
  topologyHasChanged_();
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::unRoot(bool joinRootSons)
{
  if (joinRootSons)
  {
    // the root must have exactly two joinRootSons
    std::vector<Graph::NodeId> sons = getSons(GraphImpl::getRoot());
    if (sons.size() != 2)
      throw Exception("The root must have two sons to join them.");
    GraphImpl::unlink(GraphImpl::getRoot(), sons.at(0));
    GraphImpl::unlink(GraphImpl::getRoot(), sons.at(1));
    GraphImpl::link(sons.at(0), sons.at(1));
    GraphImpl::setRoot(sons.at(0));
  }
  GraphImpl::makeUndirected();
}

template<class GraphImpl>
std::vector<Graph::NodeId> TreeGraphImpl<GraphImpl>::getSons(Graph::NodeId node) const
{
  return GraphImpl::getOutgoingNeighbors(node);
}

template<class GraphImpl>
std::vector<Graph::EdgeId> TreeGraphImpl<GraphImpl>::getBranches(Graph::NodeId node) const
{
  return GraphImpl::getOutgoingEdges(node);
}

template<class GraphImpl>
std::unique_ptr<Graph::NodeIterator> TreeGraphImpl<GraphImpl>::sonsIterator(Graph::NodeId node)
{
  return GraphImpl::outgoingNeighborNodesIterator(node);
}

template<class GraphImpl>
std::unique_ptr<Graph::EdgeIterator> TreeGraphImpl<GraphImpl>::branchesIterator(Graph::NodeId node)
{
  return GraphImpl::outgoingEdgesIterator(node);
}

template<class GraphImpl>
size_t TreeGraphImpl<GraphImpl>::getNumberOfSons(Graph::NodeId node) const
{
  return GraphImpl::getNumberOfOutgoingNeighbors(node);
}

template<class GraphImpl>
std::vector<Graph::NodeId> TreeGraphImpl<GraphImpl>::removeSons(Graph::NodeId node)
{
  std::vector<Graph::NodeId> sons = getSons(node);
  for (std::vector<Graph::NodeId>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
  {
    removeSon(node, *currSon);
  }
  return sons;
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::removeSon(Graph::NodeId node, Graph::NodeId son)
{
  GraphImpl::unlink(node, son);
  topologyHasChanged_();
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::setOutGroup(Graph::NodeId newOutGroup)
{
  mustBeRooted_();
  deleteNode(GraphImpl::getRoot());

  Graph::NodeId newRoot = GraphImpl::createNodeFromEdge(getEdge(getFather(newOutGroup), newOutGroup));
  rootAt(newRoot);
}

template<class GraphImpl>
std::vector<Graph::NodeId> TreeGraphImpl<GraphImpl>::getNodePathBetweenTwoNodes(Graph::NodeId nodeA, Graph::NodeId nodeB, bool includeAncestor) const
{
  GraphImpl::nodeMustExist_(nodeA);
  GraphImpl::nodeMustExist_(nodeB);
  std::vector<Graph::NodeId> path;
  std::vector<Graph::NodeId> pathMatrix1;
  std::vector<Graph::NodeId> pathMatrix2;

  Graph::NodeId nodeUp = nodeA;
  while (hasFather(nodeUp))
  {
    pathMatrix1.push_back(nodeUp);
    nodeUp = getFather(nodeUp);
  }
  pathMatrix1.push_back(nodeUp); // The root.

  nodeUp = nodeB;
  while (hasFather(nodeUp))
  {
    pathMatrix2.push_back(nodeUp);
    nodeUp = getFather(nodeUp);
  }
  pathMatrix2.push_back(nodeUp); // The root.
  // Must check that the two nodes have the same root!!!

  size_t tmp1 = pathMatrix1.size();
  size_t tmp2 = pathMatrix2.size();

  while ((tmp1 > 0) && (tmp2 > 0))
  {
    if (pathMatrix1[tmp1 - 1] != pathMatrix2[tmp2 - 1])
      break;
    tmp1--; tmp2--;
  }
  // (tmp1 - 1) and (tmp2 - 1) now point toward the first non-common nodes

  for (size_t y = 0; y < tmp1; ++y)
  {
    path.push_back(pathMatrix1[y]);
  }
  if (includeAncestor) // FIXME: one of the extremities may be the ancestor!!!
    path.push_back(pathMatrix1[tmp1]);                                            // pushing once, the Node that was common to both.
  for (size_t j = tmp2; j > 0; --j)
  {
    path.push_back(pathMatrix2[j - 1]);
  }
  return path;
}

template<class GraphImpl>
std::vector<Graph::EdgeId> TreeGraphImpl<GraphImpl>::getEdgePathBetweenTwoNodes(Graph::NodeId nodeA, Graph::NodeId nodeB) const
{
  std::vector<Graph::EdgeId> path;
  std::vector<Graph::NodeId> pathNodes = getNodePathBetweenTwoNodes(nodeA, nodeB, true);
  for (size_t currNodeNr = 0; currNodeNr + 1 < pathNodes.size(); currNodeNr++)
  {
    path.push_back(GraphImpl::getAnyEdge(pathNodes.at(currNodeNr), pathNodes.at(currNodeNr + 1)));
  }
  return path;
}

template<class GraphImpl>
std::vector<Graph::NodeId> TreeGraphImpl<GraphImpl>::getSubtreeNodes(Graph::NodeId localRoot) const
{
  mustBeValid_();
  std::vector<Graph::EdgeId> metNodes;
  fillSubtreeMetNodes_(metNodes, localRoot);
  return metNodes;
}

template<class GraphImpl>
std::vector<Graph::EdgeId> TreeGraphImpl<GraphImpl>::getSubtreeEdges(Graph::NodeId localRoot) const
{
  mustBeValid_();
  std::vector<Graph::EdgeId> metEdges;
  fillSubtreeMetEdges_(metEdges, localRoot);
  return metEdges;
}


template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::fillSubtreeMetNodes_(std::vector<Graph::NodeId>& metNodes, Graph::NodeId localRoot) const
{
  metNodes.push_back(localRoot);
  std::vector<Graph::NodeId> sons = GraphImpl::getOutgoingNeighbors(localRoot);
  for (std::vector<Graph::NodeId>::iterator currSon = sons.begin(); currSon != sons.end(); currSon++)
  {
    fillSubtreeMetNodes_(metNodes, *currSon);
  }
}

template<class GraphImpl>
void TreeGraphImpl<GraphImpl>::fillSubtreeMetEdges_(std::vector<Graph::EdgeId>& metEdges, Graph::NodeId localRoot) const
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
