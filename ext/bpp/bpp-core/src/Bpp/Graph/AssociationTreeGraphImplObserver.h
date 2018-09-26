//
// File AssociationTreeGraphImplObserver.h
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

#ifndef _ASSOCIATION_TREEGRAPH_IMPL_OBSERVER_HPP_
#define _ASSOCIATION_TREEGRAPH_IMPL_OBSERVER_HPP_

#include "TreeGraphImpl.h"
#include "AssociationTreeGraphObserver.h"
#include "AssociationGraphImplObserver.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <memory>

namespace bpp
{
template<class N, class E, class TreeGraphImpl>
class AssociationTreeGraphImplObserver :
  public AssociationTreeGraphObserver<N, E>,
  public AssociationGraphImplObserver<N, E, TreeGraphImpl>
{
public:
  typedef typename AssociationGraphObserver<N, E>::NodeIndex NodeIndex;
  typedef typename AssociationGraphObserver<N, E>::EdgeIndex EdgeIndex;

  typedef typename Graph::NodeId NodeGraphid;
  typedef typename Graph::EdgeId EdgeGraphid;

public:
  /**
   * Constructor
   * @param rooted if the graph rooted
   */
  AssociationTreeGraphImplObserver(bool rooted = false) :
    AssociationGraphImplObserver<N, E, TreeGraphImpl>(rooted)
  {}

  /**
   * Constructor
   * @param subjectTreeGraph the graph which is observed
   */
  AssociationTreeGraphImplObserver(std::shared_ptr<TreeGraphImpl> subjectTreeGraph = 00) :
    AssociationGraphImplObserver<N, E, TreeGraphImpl>(subjectTreeGraph)
  {}

  /**
   * Copy Constructor
   * @param treeGraphObserver the treeGraphObserver to be copied
   */

  AssociationTreeGraphImplObserver(bpp::AssociationTreeGraphImplObserver<N, E, TreeGraphImpl> const& treeGraphObserver) :
    AssociationGraphImplObserver<N, E, TreeGraphImpl>(treeGraphObserver)
  {}

  /**
   * Copy Constructor
   * @param treeGraphObserver the treeGraphObserver to be copied
   */

  template<class N2, class E2>
  AssociationTreeGraphImplObserver(bpp::AssociationTreeGraphImplObserver<N2, E2, TreeGraphImpl> const& treeGraphObserver) :
    AssociationGraphImplObserver<N, E, TreeGraphImpl>::AssociationGraphImplObserver(treeGraphObserver)
  {}


  /**
   * = Operator
   * @param treeGraphObserver the treeGraphObserver we want to copy the values
   */

  AssociationTreeGraphImplObserver<N, E, TreeGraphImpl>& operator=(bpp::AssociationTreeGraphImplObserver<N, E, TreeGraphImpl> const& treeGraphObserver)
  {
    AssociationGraphImplObserver<N, E, TreeGraphImpl>::operator=(treeGraphObserver);
    return *this;
  }


  /**
   * Destructor
   */

  ~AssociationTreeGraphImplObserver()
  {}


  /**
   * clone function
   */
  AssociationTreeGraphImplObserver<N, E, TreeGraphImpl>* clone() const { return new AssociationTreeGraphImplObserver<N, E, TreeGraphImpl>(*this); }


  /**
   * Is the graph a tree? A tree must be acyclic and with no isolated node.
   * @return true if valid tree
   */
  bool isValid() const
  {
    return this->getGraph()->isValid();
  }

  /**
   * Return, in a rooted tree, the branch leading to the father
   * @param nodeObject the concerned node
   * @return an Edge which is the branch to the father
   */
  std::shared_ptr<E>  getEdgeToFather(const std::shared_ptr<N>  nodeObject) const
  {
    return this->getEdgeFromGraphid(this->getGraph()->getEdgeToFather(this->getNodeGraphid(nodeObject)));
  }

  std::shared_ptr<E>  getEdgeToFather(const NodeIndex index) const
  {
    return this->getEdgeFromGraphid(this->getGraph()->getEdgeToFather(this->getNodeGraphid(this->getNode(index))));
  }


  /**
   * @brief Sets the root and make the tree directed from root to leaves
   *
   */
  void rootAt(const std::shared_ptr<N> root)
  {
    this->getGraph()->rootAt(this->getNodeGraphid(root));
  }

  /*
   * @brief check if rooted, ie directed
   *
   */
  bool isRooted() const
  {
    return this->getGraph()->isRooted();
  }

  /**
   * Return, in a rooted tree, the father node
   * @param nodeObject the concerned node
   * @return the father
   */

  std::shared_ptr<N>  getFather(const std::shared_ptr<N>  nodeObject) const
  {
    return this->getNodeFromGraphid(this->getGraph()->getFather(this->getNodeGraphid(nodeObject)));
  }

  /**
   * Has the node a father?
   */
  bool hasFather(const std::shared_ptr<N>  nodeObject) const
  {
    return this->getGraph()->hasFather(this->getNodeGraphid(nodeObject));
  }

  bool hasFather(const NodeIndex index) const
  {
    return this->getGraph()->hasFather(this->getNodeGraphid(this->getNode(index)));
  }

  /**
   * Return, in a rooted tree, the sons of a node
   * @param node the concerned node
   * @return a vector of son Nodes
   */

  std::vector<std::shared_ptr<N> > getSons(const std::shared_ptr<N>  node) const
  {
    return this->getNodesFromGraphid(this->getGraph()->getSons(this->getNodeGraphid(node)));
  }

  std::vector<NodeIndex> getSons(const NodeIndex node) const
  {
    return this->getNodeIndexes(this->getNodesFromGraphid(this->getGraph()->getSons(this->getNodeGraphid(this->getNode(node)))));
  }

  /**
   * Return, in a rooted tree, the branches to the sons of a node
   * @param node the concerned node
   * @return a vector of branch Nodes
   */

  std::vector<std::shared_ptr<E> > getBranches(const std::shared_ptr<N>  node) const
  {
    return this->getEdgesFromGraphid(this->getGraph()->getBranches(this->getNodeGraphid(node)));
  }

  std::vector<EdgeIndex> getBranches(const NodeIndex node) const
  {
    return this->getEdgeIndexes(this->getEdgesFromGraphid(this->getGraph()->getBranches(this->getNodeGraphid(this->getNode(node)))));
  }

  /**
   * Return, in a rooted tree, the son of an edge
   * @param edge the concerned edge
   * @return the son Node
   */

  std::shared_ptr<N> getSon(const std::shared_ptr<E>  edge) const
  {
    return this->getNodeFromGraphid(this->getGraph()->getBottom(this->getEdgeGraphid(edge)));
  }

  NodeIndex getSon(const EdgeIndex edge) const
  {
    return this->getNodeIndex(this->getNode(this->getGraph()->getBottom(this->getEdgeGraphid(this->getEdge(edge)))));
  }

  /**
   * Return, in a rooted tree, the father of an edge
   * @param edge the concerned edge
   * @return the father Node
   */

  std::shared_ptr<N> getFather(const std::shared_ptr<E>  edge) const
  {
    return this->getNodeFromGraphid(this->getGraph()->getTop(this->getEdgeGraphid(edge)));
  }

  NodeIndex getFather(const EdgeIndex edge) const
  {
    return this->getNodeIndex(this->getNode(this->getGraph()->getTop(this->getEdgeGraphid(this->getEdge(edge)))));
  }


  /**
   * Return, in a rooted tree, the number of sons
   * @param node the concerned node
   * @return the number of sons
   */
  size_t getNumberOfSons(const std::shared_ptr<N>  node) const
  {
    return this->getGraph()->getNumberOfSons(this->getNodeGraphid(node));
  }

  /**
   * Get the leaves of a tree under a peculiar node.
   *
   * @param node the starting node
   * @return a vector containing the leaves
   */

  std::vector<std::shared_ptr<N> > getLeavesUnderNode(std::shared_ptr<N>  node) const
  {
    return this->getNodesFromGraphid(this->getGraph()->getLeavesUnderNode(this->getNodeGraphid(node)));
  }


  /**
   * Remove the sons of a node
   * @return a vector containing the removed nodes
   */

  std::vector<std::shared_ptr<N> > removeSons(const std::shared_ptr<N>  node)
  {
    return this->getNodesFromGraphid(this->getGraph()->removeSons(this->getNodeGraphid(node)));
  }

  /**
   * Remove a son of a node
   */
  void removeSon(const std::shared_ptr<N> node, const std::shared_ptr<N> son)
  {
    this->getGraph()->removeSon(this->getNodeGraphid(node), this->getNodeGraphid(son));
  }

  /**
   * Change / set the father of a node
   * @param nodeObject the concerned node
   * @param fatherNodeObject the node to be the father
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */
  void setFather(const std::shared_ptr<N>  nodeObject, const std::shared_ptr<N> fatherNodeObject, const std::shared_ptr<E> edgeObject = 0)
  {
    if (edgeObject)
      this->getGraph()->setFather(this->getNodeGraphid(nodeObject), this->getNodeGraphid(fatherNodeObject), this->getEdgeGraphid(edgeObject));
    else
      this->getGraph()->setFather(this->getNodeGraphid(nodeObject), this->getNodeGraphid(fatherNodeObject));
  }


  /**
   * Add a son to a node
   * @param nodeObject the concerned node
   * @param sonNodeObject the node to be added as a son to the father
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */
  void addSon(const std::shared_ptr<N>  nodeObject, const std::shared_ptr<N> sonNodeObject, const std::shared_ptr<E> edgeObject = 0)
  {
    if (edgeObject)
      this->getGraph()->addSon(this->getNodeGraphid(nodeObject), this->getNodeGraphid(sonNodeObject), this->getEdgeGraphid(edgeObject));
    else
      this->getGraph()->addSon(this->getNodeGraphid(nodeObject), this->getNodeGraphid(sonNodeObject));
  }

  /**
   * Iterators
   *
   */

  /*
   * @brief builds iterator on the sons of a Node
   *
   */

  std::unique_ptr<typename AssociationTreeGraphObserver<N, E>::NodeIterator> sonsIterator(std::shared_ptr<N> node)
  {
    return this->outgoingNeighborNodesIterator(node);
  }

  /*
   * @brief builds iterator on the branches to sons of a Node
   *
   */

  std::unique_ptr<typename AssociationTreeGraphObserver<N, E>::EdgeIterator> branchesIterator(std::shared_ptr<N> node)
  {
    return this->outgoingEdgesIterator(node);
  }

  /**
   * @brief Get a vector of ancestor nodes between to nodes.
   *
   * @param nodeA first node.
   * @param nodeB second node.
   * @param includeAncestor Tell if the common ancestor must be included in the vector.
   * @return A vector of ancestor nodes ids.
   * @throw PhyloNodeNotFoundException If a node is not found.
   */

  std::vector<std::shared_ptr<N> > getNodePathBetweenTwoNodes(const std::shared_ptr<N>  nodeA, const std::shared_ptr<N>  nodeB, bool includeAncestor = true) const
  {
    return this->getNodesFromGraphid(this->getGraph()->getNodePathBetweenTwoNodes(this->getNodeGraphid(nodeA), this->getNodeGraphid(nodeB), includeAncestor));
  }

  std::vector<std::shared_ptr<E> > getEdgePathBetweenTwoNodes(const std::shared_ptr<N>  nodeA, const std::shared_ptr<N>  nodeB) const
  {
    return this->getEdgesFromGraphid(this->getGraph()->getEdgePathBetweenTwoNodes(this->getNodeGraphid(nodeA), this->getNodeGraphid(nodeB)));
  }

  std::vector<std::shared_ptr<N> > getSubtreeNodes(const std::shared_ptr<N> localRoot)
  {
    return this->getNodesFromGraphid(this->getGraph()->getSubtreeNodes(this->getNodeGraphid(localRoot)));
  }

  std::vector<std::shared_ptr<E> > getSubtreeEdges(const std::shared_ptr<N> localRoot)
  {
    return AssociationGraphImplObserver<N, E, TreeGraphImpl>::getEdgesFromGraphid(this->getGraph()->getSubtreeEdges(this->getNodeGraphid(localRoot)));
  }
};

/********************/

template<class N, class E>
using AssociationTreeGlobalGraphObserver =  AssociationTreeGraphImplObserver<N, E, TreeGlobalGraph>;
}


#endif
