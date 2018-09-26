//
// File AssociationTreeGraphObserver.h
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

#ifndef _ASSOCIATION_TREE_GRAPHOBSERVER_HPP_
#define _ASSOCIATION_TREE_GRAPHOBSERVER_HPP_

#include "TreeGraph.h"
#include "AssociationGraphObserver.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <memory>

namespace bpp
{
/**
 * @brief Defines a Tree Graph Associator. It is a template which follows
 * (subscribed to) a Graph.
 *
 * @author Thomas Bigot
 */

template<class N, class E>
class AssociationTreeGraphObserver :
  public virtual AssociationGraphObserver<N, E>
{
public:
  typedef typename AssociationTreeGraphObserver<N, E>::NodeIndex NodeIndex;
  typedef typename AssociationTreeGraphObserver<N, E>::EdgeIndex EdgeIndex;

  typedef typename TreeGraph::NodeId NodeGraphid;
  typedef typename TreeGraph::EdgeId EdgeGraphid;

public:
  /**
   * Is the graph a tree? A tree must be acyclic and with no isolated node.
   * @return true if valid tree
   */

  virtual bool isValid() const = 0;

  /**
   * Return, in a rooted tree, the branch leading to the father
   * @param nodeObject the concerned node
   * @return an Edge which is the branch to the father
   */

  virtual std::shared_ptr<E>  getEdgeToFather(const std::shared_ptr<N>  nodeObject) const = 0;

  virtual std::shared_ptr<E>  getEdgeToFather(const NodeIndex nodeIndex) const = 0;

  /**
   * @brief Sets the root and make the tree directed from root to leaves
   *
   */

  virtual void rootAt(const std::shared_ptr<N> root) = 0;

  /*
   * @brief check if rooted, ie directed
   *
   */

  virtual bool isRooted() const = 0;

  /**
   * Return, in a rooted tree, the father node
   * @param nodeObject the concerned node
   * @return the father
   */

  virtual std::shared_ptr<N>  getFather(const std::shared_ptr<N>  nodeObject) const = 0;

  /**
   * Has the node a father?
   */

  virtual bool hasFather(const std::shared_ptr<N>  nodeObject) const = 0;
  virtual bool hasFather(const NodeIndex node) const = 0;

  /**
   * Return, in a rooted tree, the sons of a node
   * @param node the concerned node
   * @return a vector of son Nodes
   */

  virtual std::vector<std::shared_ptr<N> > getSons(const std::shared_ptr<N>  node) const = 0;

  virtual std::vector<NodeIndex> getSons(const NodeIndex node) const = 0;

  /**
   * Return the branches to the sons of a node
   * @param node the concerned node
   * @return a vector of branch Nodes
   */

  virtual std::vector<std::shared_ptr<E> > getBranches(const std::shared_ptr<N>  node) const = 0;

  virtual std::vector<EdgeIndex> getBranches(const NodeIndex node) const = 0;

  /**
   * Return, in a rooted tree, the son of an edge
   * @param edge the concerned edge
   * @return the son Node
   */

  virtual std::shared_ptr<N> getSon(const std::shared_ptr<E>  edge) const = 0;
  virtual NodeIndex getSon(const EdgeIndex edge) const = 0;

  /**
   * Return, in a rooted tree, the father of an edge
   * @param edge the concerned edge
   * @return the father Node
   */

  virtual std::shared_ptr<N> getFather(const std::shared_ptr<E>  edge) const = 0;
  virtual NodeIndex getFather(const EdgeIndex edge) const = 0;

  /**
   * Return, in a rooted tree, the number of sons
   * @param node the concerned node
   * @return the number of sons
   */

  virtual size_t getNumberOfSons(const std::shared_ptr<N>  node) const = 0;

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
   */

  virtual void removeSon(const std::shared_ptr<N> node, const std::shared_ptr<N> son) = 0;

  /**
   * Change / set the father of a node
   * @param nodeObject the concerned node
   * @param fatherNodeObject the node to be the father
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */

  virtual void setFather(const std::shared_ptr<N>  nodeObject, const std::shared_ptr<N>  fatherNodeObject, const std::shared_ptr<E> edgeObject = 0) = 0;

  /**
   * Add a son to a node
   * @param nodeObject the concerned node
   * @param sonNodeObject the node to be added as a son to the father
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */

  virtual void addSon(const std::shared_ptr<N>  nodeObject, const std::shared_ptr<N>  sonNodeObject, const std::shared_ptr<E> edgeObject = 0) = 0;

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
   * @brief builds iterator on the branches to sons of a Node
   *
   */

  virtual std::unique_ptr<EdgeIterator> branchesIterator(std::shared_ptr<N> node) = 0;

  /**
   * @brief Get a vector of ancestor nodes between to nodes.
   *
   * @param nodeA first node.
   * @param nodeB second node.
   * @param includeAncestor Tell if the common ancestor must be included in the vector.
   * @return A vector of ancestor nodes ids.
   * @throw PhyloNodeNotFoundException If a node is not found.
   */

  virtual std::vector<std::shared_ptr<N> > getNodePathBetweenTwoNodes(const std::shared_ptr<N>  nodeA, const std::shared_ptr<N>  nodeB, bool includeAncestor = true) const = 0;

  virtual std::vector<std::shared_ptr<E> > getEdgePathBetweenTwoNodes(const std::shared_ptr<N>  nodeA, const std::shared_ptr<N>  nodeB) const = 0;

  virtual std::vector<std::shared_ptr<N> > getSubtreeNodes(const std::shared_ptr<N> localRoot) = 0;

  virtual std::vector<std::shared_ptr<E> > getSubtreeEdges(const std::shared_ptr<N> localRoot) = 0;
};
}

#endif
