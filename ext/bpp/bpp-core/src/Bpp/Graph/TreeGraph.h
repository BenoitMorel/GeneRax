//
// File TreeGraph.h
// Created by: Thomas Bigot
// Last modification : vendredi 4 novembre 2016, à 10h 25
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


#ifndef _TREEGRAPH_H_
#define _TREEGRAPH_H_

#include <string>
#include <vector>
#include <iostream>
#include <ostream>


#include "Graph.h"

#include "../Exceptions.h"
#include "../Numeric/VectorTools.h"

namespace bpp
{
class TreeGraph :
  public virtual Graph
{
public:
  /**
   * Is the graph a tree? A tree must be acyclic and with no isolated node.
   * @return true if valid tree
   */
  virtual bool isValid() const = 0;

  /**
   * Get the father node of a node in a rooted tree
   * @return the father node
   */

  virtual Graph::NodeId getFather(Graph::NodeId nodeid) const = 0;

  /**
   * Get the branch leading to the father in a rooted tree
   * @return the branch between a node and its father
   */

  virtual Graph::EdgeId getEdgeToFather(Graph::NodeId node) const = 0;

  /**
   * Check if node has a father
   */

  virtual bool hasFather(Graph::NodeId node) const = 0;

  /**
   * Says if  a node is a leaf (ie has at most one neighbor).
   */

  virtual bool isLeaf(Graph::NodeId node) const = 0;

  /**
   * Get the leaves under a node
   * @param node the starting node
   * @return a vector containing the leaves
   */

  virtual std::vector<Graph::NodeId> getLeavesUnderNode(Graph::NodeId node) const = 0;

  /**
   * Get the sons node of a node
   */

  virtual std::vector<Graph::NodeId> getSons(Graph::NodeId node) const = 0;

  /**
   * Get the branches to the sons of a node
   */

  virtual std::vector<Graph::EdgeId> getBranches(Graph::NodeId node) const = 0;

  /**
   * Get a iterator on the sons node of a node
   */

  virtual std::unique_ptr<Graph::NodeIterator> sonsIterator(Graph::NodeId node) = 0;

  /**
   * Get a iterator on the branches to sons of a node
   */

  virtual std::unique_ptr<Graph::EdgeIterator> branchesIterator(Graph::NodeId node) = 0;

  /**
   * @brief Get the number of sons node
   */

  virtual size_t getNumberOfSons(Graph::NodeId node) const = 0;

  /**
   * set the father node of a node in a rooted tree
   */

  virtual void setFather(Graph::NodeId node, Graph::NodeId fatherNode) = 0;

  virtual void setFather(Graph::NodeId node, Graph::NodeId fatherNode, Graph::EdgeId edgeId) = 0;

  /**
   * Add a son to a node in a rooted tree
   */

  virtual void addSon(Graph::NodeId node, Graph::NodeId sonNode) = 0;

  virtual void addSon(Graph::NodeId node, Graph::NodeId sonNode, Graph::EdgeId edgeId) = 0;

  /**
   * Remove all the sons
   */

  std::vector<Graph::NodeId> removeSons(Graph::NodeId node);

  /**
   * Remove one son
   */

  virtual void removeSon(Graph::NodeId node, Graph::NodeId son) = 0;

  /**
   * Re-root the tree with the new root
   */

  virtual void rootAt(Graph::NodeId newRoot) = 0;

  /**
   * Set the tree to its flat unrooted version.
   * As an algorithmical convenience, a root node is kept, but it has
   * no logical significance.
   */

  virtual void unRoot(bool joinRootSons) = 0;

  /**
   * Set a node as a new outgroup in a rooted tree, will make a root between
   * the given node and its father.
   */

  virtual void setOutGroup(Graph::NodeId newOutGroup) = 0;

  /**
   * Get all the nodes of a subtree
   */

  virtual std::vector<Graph::NodeId> getSubtreeNodes(Graph::NodeId localRoot) const = 0;

  /**
   * Get all the branches of a subtree
   */

  virtual std::vector<Graph::EdgeId> getSubtreeEdges(Graph::NodeId localRoot) const = 0;

  // ///FROM TREETOOLS & TREETOOLS COMPAT


  virtual std::vector<Graph::NodeId> getNodePathBetweenTwoNodes(Graph::NodeId nodeA, Graph::NodeId nodeB, bool includeAncestor = true) const = 0;
  virtual std::vector<Graph::EdgeId> getEdgePathBetweenTwoNodes(Graph::NodeId nodeA, Graph::NodeId nodeB) const = 0;
};
}


#endif
