//
// File AssociationDAGraphObserver.h
// Created by: Laurent Guéguen
// Last modification : lundi 19 décembre 2016, à 21h 48
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

#ifndef _ASSOCIATION_DA_GRAPHOBSERVER_HPP_
#define _ASSOCIATION_DA_GRAPHOBSERVER_HPP_

#include "DAGraph.h"
#include "AssociationGraphObserver.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <memory>

namespace bpp
{
/**
 * @brief Defines a DA Graph Associator. It is a template which
 * follows (subscribed to) a Graph.
 *
 * @author Laurent Guéguen
 */

template<class N, class E>
class AssociationDAGraphObserver :
  public virtual AssociationGraphObserver<N, E>
{
public:
  typedef typename AssociationDAGraphObserver<N, E>::NodeIndex NodeIndex;
  typedef typename AssociationDAGraphObserver<N, E>::EdgeIndex EdgeIndex;

  typedef typename DAGraph::NodeId NodeGraphid;
  typedef typename DAGraph::EdgeId EdgeGraphid;

public:
  /**
   * Is the graph a DAG? A DAG must be acyclic.
   * @return true if valid DAG
   */

  virtual bool isValid() const = 0;

  /**
   * Is the DAG rooted? Ie with at most one node with no father.
   * @return true if rooted DAG
   */

  virtual bool isRooted() const = 0;

  /**
   * @brief Sets the root and make the DAG directed from root to
   * leaves
   *
   */

  virtual void rootAt(const std::shared_ptr<N> root) = 0;


  /**
   * Return, the fathers of a node
   * @param node the concerned node
   * @return the father
   */

  virtual std::vector<std::shared_ptr<N> > getFathers(const std::shared_ptr<N>  node) const = 0;

  virtual std::vector<NodeIndex> getFathers(const NodeIndex node) const = 0;

  /**
   * Has the node a father?
   */

  virtual bool hasFather(const std::shared_ptr<N>  nodeObject) const = 0;
  virtual bool hasFather(const NodeIndex node) const = 0;

  /**
   * Return the sons of a node
   * @param node the concerned node
   * @return a vector of son Nodes
   */

  virtual std::vector<std::shared_ptr<N> > getSons(const std::shared_ptr<N>  node) const = 0;

  virtual std::vector<NodeIndex> getSons(const NodeIndex node) const = 0;

  /**
   * Return the son of an edge
   * @param edge the concerned edge
   * @return the son Node
   */

  virtual std::shared_ptr<N> getSon(const std::shared_ptr<E>  edge) const = 0;
  virtual NodeIndex getSon(const EdgeIndex edge) const = 0;

  /**
   * Return the father of an edge
   * @param edge the concerned edge
   * @return the father Node
   */

  virtual std::shared_ptr<N> getFather(const std::shared_ptr<E>  edge) const = 0;
  virtual NodeIndex getFather(const EdgeIndex edge) const = 0;

  /**
   * Return the number of sons
   * @param node the concerned node
   * @return the number of sons
   */

  virtual size_t getNumberOfSons(const std::shared_ptr<N>  node) const = 0;

  /**
   * Return the number of fathers
   * @param node the concerned node
   * @return the number of fathers
   */

  virtual size_t getNumberOfFathers(const std::shared_ptr<N>  node) const = 0;

  /**
   * Get the leaves of a tree under a peculiar node.
   *
   * @param node the starting node
   * @return a vector containing the leaves
   */

  virtual std::vector<std::shared_ptr<N> > getLeavesUnderNode(std::shared_ptr<N>  node) const = 0;

  /**
   * Remove the sons of a node
   * @return a vector containing the removed nodes
   */

  virtual std::vector<std::shared_ptr<N> > removeSons(const std::shared_ptr<N>  node) = 0;

  /**
   * Remove a son of a node
   * @param node the concerned node
   * @param son the node to be removed
   */

  virtual void removeSon(const std::shared_ptr<N> node, const std::shared_ptr<N> son) = 0;

  /**
   * Add a son to a node
   * @param nodeObject the concerned node
   * @param sonNodeObject the node to be added as a son to the father
   * @param edgeObject the optional edge  between the nodes (default
   * = 0)
   */

  virtual void addSon(const std::shared_ptr<N>  nodeObject, const std::shared_ptr<N>  sonNodeObject, const std::shared_ptr<E> edgeObject = 0) = 0;

  /**
   * Remove the fathers of a node
   * @return a vector containing the removed fathers
   */

  virtual std::vector<std::shared_ptr<N> > removeFathers(const std::shared_ptr<N> node) = 0;

  /**
   * Remove a father of a node
   */

  virtual void removeFather(const std::shared_ptr<N> node, const std::shared_ptr<N> father) = 0;

  /**
   * Add a father to a node
   * @param nodeObject the concerned node
   * @param fatherNodeObject the node to be added as a father to the node
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */

  virtual void addFather(const std::shared_ptr<N>  nodeObject, const std::shared_ptr<N>  fatherNodeObject, const std::shared_ptr<E> edgeObject = 0) = 0;

  /**
   * Iterators
   *
   */

  /*
   * @brief builds iterator on the neighbor nodes of a Node
   *
   */

  typedef typename AssociationGraphObserver<N, E>::NodeIterator NodeIterator;
  typedef typename AssociationGraphObserver<N, E>::EdgeIterator EdgeIterator;


  /*
   * @brief builds iterator on the sons of a Node
   *
   */

  virtual std::unique_ptr<NodeIterator> sonsIterator(std::shared_ptr<N> node) = 0;

  /*
   * @brief builds iterator on the fathers of a Node
   *
   */

  virtual std::unique_ptr<NodeIterator> fathersIterator(std::shared_ptr<N> node) = 0;

  /**
   * @brief Get Below Objects.
   *
   * @param localRoot of the subgraph.
   * @return A vector of ancestor nodes ids.
   *
   */

  virtual std::vector<std::shared_ptr<N> > getBelowNodes(const std::shared_ptr<N> localRoot) = 0;

  virtual std::vector<std::shared_ptr<E> > getBelowEdges(const std::shared_ptr<N> localRoot) = 0;
};
}

#endif
