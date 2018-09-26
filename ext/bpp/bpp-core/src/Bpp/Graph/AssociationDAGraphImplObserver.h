//
// File AssociationDAGraphImplObserver.h
// Created by: Laurent Guéguen
// Last modification : lundi 19 décembre 2016, à 22h 14
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

#ifndef _ASSOCIATION_DAGRAPH_IMPL_OBSERVER_HPP_
#define _ASSOCIATION_DAGRAPH_IMPL_OBSERVER_HPP_

#include "DAGraphImpl.h"
#include "AssociationDAGraphObserver.h"
#include "AssociationGraphImplObserver.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>
#include <memory>

namespace bpp
{
template<class N, class E, class DAGraphImpl>
class AssociationDAGraphImplObserver :
  public AssociationDAGraphObserver<N, E>,
  public AssociationGraphImplObserver<N, E, DAGraphImpl>
{
public:
  typedef typename AssociationGraphObserver<N, E>::NodeIndex NodeIndex;
  typedef typename AssociationGraphObserver<N, E>::EdgeIndex EdgeIndex;

  typedef typename Graph::NodeId NodeGraphid;
  typedef typename Graph::EdgeId EdgeGraphid;

public:
  /**
   * Constructor
   */

  AssociationDAGraphImplObserver() :
    AssociationGraphImplObserver<N, E, DAGraphImpl>(true)
  {}

  /**
   * Constructor
   * @param subjectDAGraph the DAGraph which is observed
   */

  AssociationDAGraphImplObserver(std::shared_ptr<DAGraphImpl> subjectDAGraph) :
    AssociationGraphImplObserver<N, E, DAGraphImpl>(subjectDAGraph)
  {}

  /**
   * Copy Constructor
   * @param dAGraphObserver the DAGraphObserver to be copied
   */

  AssociationDAGraphImplObserver(bpp::AssociationDAGraphImplObserver<N, E, DAGraphImpl> const& dAGraphObserver) :
    AssociationGraphImplObserver<N, E, DAGraphImpl>(dAGraphObserver)
  {}


  /**
   * Copy Constructor
   * @param dAGraphObserver the DAGraphObserver to be copied
   */

  template<class N2, class E2>
  AssociationDAGraphImplObserver(bpp::AssociationDAGraphImplObserver<N2, E2, DAGraphImpl> const& dAGraphObserver) :
    AssociationGraphImplObserver<N, E, DAGraphImpl>(dAGraphObserver)
  {}

  /**
   * = Operator
   * @param dAGraphObserver the DAGraphObserver we want to copy the values
   */

  AssociationDAGraphImplObserver<N, E, DAGraphImpl>& operator=(bpp::AssociationDAGraphImplObserver<N, E, DAGraphImpl> const& dAGraphObserver)
  {
    AssociationGraphImplObserver<N, E, DAGraphImpl>::operator=(dAGraphObserver);
    return *this;
  }

  /**
   * Destructor
   */

  ~AssociationDAGraphImplObserver()
  {}

  /**
   * clone function
   */

  AssociationDAGraphImplObserver<N, E, DAGraphImpl>* clone() const
  {
    return new AssociationDAGraphImplObserver<N, E, DAGraphImpl>(*this);
  }


  /**
   * Is the graph a DAG? A DAG must be acyclic and with no isolated node.
   * @return true if valid DAG
   */
  bool isValid() const
  {
    return this->getGraph()->isValid();
  }

  /**
   * Is the DAG rooted? Ie with a unique node with no father.
   * @return true if rooted
   */
  bool isRooted() const
  {
    return this->getGraph()->isRooted();
  }


  /**
   * Return the fathers of a node
   * @param node the concerned node
   * @return a vector of son Nodes
   */

  std::vector<std::shared_ptr<N> > getFathers(const std::shared_ptr<N>  node) const
  {
    return this->getNodesFromGraphid(this->getGraph()->getFathers(this->getNodeGraphid(node)));
  }

  std::vector<NodeIndex> getFathers(const NodeIndex node) const
  {
    return this->getNodeIndexes(this->getNodesFromGraphid(this->getGraph()->getFathers(this->getNodeGraphid(this->getNode(node)))));
  }

  /**
   * @brief Sets the root and make the DAG directed from root to
   * leaves
   *
   */
  void rootAt(const std::shared_ptr<N> root)
  {
    this->getGraph()->rootAt(this->getNodeGraphid(root));
  }

  /**
   * Has the node a father?
   */
  bool hasFather(const std::shared_ptr<N> nodeObject) const
  {
    return this->getGraph()->hasFather(this->getNodeGraphid(nodeObject));
  }

  bool hasFather(const NodeIndex node) const
  {
    return this->getGraph()->hasFather(this->getNodeGraphid(this->getNode(node)));
  }

  /**
   * Return the sons of a node
   * @return a vector of son Nodes
   */

  std::vector<std::shared_ptr<N> > getSons(const std::shared_ptr<N> node) const
  {
    return this->getNodesFromGraphid(this->getGraph()->getSons(this->getNodeGraphid(node)));
  }

  std::vector<NodeIndex> getSons(const NodeIndex node) const
  {
    return this->getNodeIndexes(this->getNodesFromGraphid(this->getGraph()->getSons(this->getNodeGraphid(this->getNode(node)))));
  }


  /**
   * Return the son of an edge
   * @param edge the concerned edge
   * @return the son Node
   */

  std::shared_ptr<N> getSon(const std::shared_ptr<E>  edge) const
  {
    return this->getNodeFromGraphid(this->getGraph()->getBottom(this->getEdgeGraphid(edge)));
  }

  NodeIndex getSon(const EdgeIndex edgeId) const
  {
    return this->getNodeIndex(this->getNode(this->getGraph()->getBottom(this->getEdgeGraphid(this->getEdge(edgeId)))));
  }


  /**
   * Return the father of an edge
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
   * Return the number of fathers
   * @param node the concerned node
   * @return the number of fathers
   */
  size_t getNumberOfFathers(const std::shared_ptr<N>  node) const
  {
    return this->getGraph()->getNumberOfFathers(this->getNodeGraphid(node));
  }

  /**
   * Return the number of sons
   * @param node the concerned node
   * @return the number of sons
   */
  size_t getNumberOfSons(const std::shared_ptr<N>  node) const
  {
    return this->getGraph()->getNumberOfSons(this->getNodeGraphid(node));
  }

  /**
   * Get the leaves of a DAG under a peculiar node.
   *
   * @param node the starting node
   * @return a vector containing the leaves
   */

  std::vector<std::shared_ptr<N> > getLeavesUnderNode(std::shared_ptr<N>  node) const
  {
    return this->getNodesFromGraphid(this->getGraph()->getLeavesUnderNode(this->getNodeGraphid(node)));
  }

  /**
   * Remove the fathers of a node
   * @return a vector containing the removed fathers
   */

  std::vector<std::shared_ptr<N> > removeFathers(const std::shared_ptr<N>  node)
  {
    return this->getNodesFromGraphid(this->getGraph()->removeFathers(this->getNodeGraphid(node)));
  }

  /**
   * Remove a father of a node
   * @return a vector containing the removed nodes
   */
  void removeFather(const std::shared_ptr<N> node, const std::shared_ptr<N> father)
  {
    this->getGraph()->removeFather(this->getNodeGraphid(node), this->getNodeGraphid(father));
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
   * Add a father to a node
   * @param nodeObject the concerned node
   * @param fatherNodeObject the node to be added as a father to the node
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */
  void addFather(const std::shared_ptr<N> nodeObject, const std::shared_ptr<N>  fatherNodeObject, const std::shared_ptr<E> edgeObject = 0)
  {
    if (edgeObject)
      try
      {
        this->getGraph()->addFather(this->getNodeGraphid(nodeObject), this->getNodeGraphid(fatherNodeObject), this->getEdgeGraphid(edgeObject));
      }
      catch (Exception& e)
      {
        this->link(fatherNodeObject, nodeObject, edgeObject);
      }
    else
      this->getGraph()->addFather(this->getNodeGraphid(nodeObject), this->getNodeGraphid(fatherNodeObject));
  }

  /**
   * Add a son to a node
   * @param nodeObject the concerned node
   * @param sonNodeObject the node to be added as a son to the father
   * @param edgeObject the optional edge  between the nodes (default
   * = 00)
   */
  void addSon(const std::shared_ptr<N> nodeObject, const std::shared_ptr<N> sonNodeObject, const std::shared_ptr<E> edgeObject = 0)
  {
    if (edgeObject)
      try
      {
        this->getGraph()->addSon(this->getNodeGraphid(nodeObject), this->getNodeGraphid(sonNodeObject), this->getEdgeGraphid(edgeObject));
      }
      catch (Exception& e)
      {
        this->link(nodeObject, sonNodeObject, edgeObject);
      }
    else
      this->getGraph()->addSon(this->getNodeGraphid(nodeObject), this->getNodeGraphid(sonNodeObject));
  }

  /**
   * Iterators
   *
   */

  /*
   * @brief builds iterator on the fathers of a Node
   *
   */

  std::unique_ptr<typename AssociationDAGraphObserver<N, E>::NodeIterator> fathersIterator(std::shared_ptr<N> node)
  {
    return this->incomingNeighborNodesIterator(node);
  }

  /*
   * @brief builds iterator on the sons of a Node
   *
   */

  std::unique_ptr<typename AssociationDAGraphObserver<N, E>::NodeIterator> sonsIterator(std::shared_ptr<N> node)
  {
    return this->outgoingNeighborNodesIterator(node);
  }

  /**
   * @brief Get Below Objects.
   *
   * @param localRoot of the subgraph.
   * @return A vector of ancestor nodes ids.
   *
   */

  std::vector<std::shared_ptr<N> > getBelowNodes(const std::shared_ptr<N> localRoot)
  {
    return this->getNodesFromGraphid(this->getGraph()->getBelowNodes(this->getNodeGraphid(localRoot)));
  }

  std::vector<std::shared_ptr<E> > getBelowEdges(const std::shared_ptr<N> localRoot)
  {
    return AssociationGraphImplObserver<N, E, DAGraphImpl>::getEdgesFromGraphid(this->getGraph()->getBelowEdges(this->getNodeGraphid(localRoot)));
  }
};

/********************/

template<class N, class E>
using AssociationDAGlobalGraphObserver =  AssociationDAGraphImplObserver<N, E, DAGlobalGraph>;

/********************/
}


#endif
