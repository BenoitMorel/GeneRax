//
// File DAGraph.h
// Created by: Laurent  Guéguen
// Last modification : lundi 19 décembre 2016, à 22h 46
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

#ifndef _DA_GRAPH_H_
#define _DA_GRAPH_H_

#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <ostream>


#include "Graph.h"

#include "../Exceptions.h"
#include "../Numeric/VectorTools.h"

namespace bpp
{
class DAGraph :
  public virtual Graph
{
public:
  /**
   * Is the graph a directed acyclic?
   * @return true if valid DAG
   */
  virtual bool isValid() const = 0;

  /**
   * Is the DAG rooted?
   *
   * @return true if rooted, ie has only one node with no father.
   */

  virtual bool isRooted() const = 0;

  /**
   * Check if node has a father
   */

  virtual bool hasFather(Graph::NodeId node) const = 0;

  /**
   * Get the fathers node of a node
   */

  virtual std::vector<Graph::NodeId> getFathers(Graph::NodeId node) const = 0;

  /**
   * @brief Get the number of fathers nodes
   */

  virtual size_t getNumberOfFathers(Graph::NodeId node) const = 0;

  /**
   * Add a father to a node
   */

  virtual void addFather(Graph::NodeId node, Graph::NodeId father) = 0;

  virtual void addFather(Graph::NodeId node, Graph::NodeId father, Graph::EdgeId edgeId) = 0;

  /**
   * Remove one father
   */

  virtual void removeFather(Graph::NodeId node, Graph::NodeId father) = 0;

  /**
   * Remove all the fathers
   */

  virtual std::vector<Graph::NodeId> removeFathers(Graph::NodeId node) = 0;

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
   * @brief Get the number of sons node
   */

  virtual size_t getNumberOfSons(Graph::NodeId node) const = 0;

  /**
   * Add a son to a node
   */

  void addSon(Graph::NodeId node, Graph::NodeId sonNode);

  void addSon(Graph::NodeId node, Graph::NodeId sonNode, Graph::EdgeId edge);

  /**
   * Remove all the sons
   */

  virtual std::vector<Graph::NodeId> removeSons(Graph::NodeId node) = 0;

  /**
   * Remove one son
   */

  virtual void removeSon(Graph::NodeId node, Graph::NodeId son) = 0;

  /**
   * Re-root the DA with the new root
   */

  virtual void rootAt(Graph::NodeId newRoot) = 0;

  /**
   * Get all the nodes below a node
   */

  virtual std::vector<Graph::NodeId> getBelowNodes(Graph::NodeId localRoot) const = 0;

  /**
   * Get all the branches below a node
   */

  virtual std::vector<Graph::EdgeId> getBelowEdges(Graph::NodeId localRoot) const = 0;
};
}


#endif
