//
// File GlobalGraph.cpp
// Created by: Thomas Bigot
// Last modification : vendredi 4 novembre 2016, à 10h 22
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


#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include "GraphObserver.h"
#include "GlobalGraph.h"
#include "../Exceptions.h"

using namespace bpp;
using namespace std;

GlobalGraph::GlobalGraph(bool directed_p) :
  directed_(directed_p),
  observers_(set<GraphObserver*>()),
  highestNodeID_(0),
  highestEdgeID_(0),
  nodeStructure_(nodeStructureType()),
  edgeStructure_(edgeStructureType()),
  root_(0)
{}


GlobalGraph::GlobalGraph(const GlobalGraph& gg) :
  directed_(gg.directed_),
  observers_(gg.observers_),
  highestNodeID_(gg.highestNodeID_),
  highestEdgeID_(gg.highestEdgeID_),
  nodeStructure_(gg.nodeStructure_),
  edgeStructure_(gg.edgeStructure_),
  root_(gg.root_)
{}

GlobalGraph& GlobalGraph::operator=(const GlobalGraph& gg)
{
  directed_ = gg.directed_;
  observers_ = gg.observers_;
  highestNodeID_ = gg.highestNodeID_;
  highestEdgeID_ = gg.highestEdgeID_;
  nodeStructure_ = gg.nodeStructure_;
  edgeStructure_ = gg.edgeStructure_;
  root_ = gg.root_;

  return *this;
}


void GlobalGraph::nodeMustExist_(const GlobalGraph::Node& node, string name) const
{
  if (nodeStructure_.find(node) == nodeStructure_.end())
  {
    ostringstream errMessage;
    errMessage << "This node must exist: " << node << " as " << name << ".";
    throw (Exception(errMessage.str()));
  }
}

void GlobalGraph::edgeMustExist_(const GlobalGraph::Edge& edge, string name) const
{
  if (edgeStructure_.find(edge) == edgeStructure_.end())
  {
    ostringstream errMessage;
    errMessage << "This edge must exist: " << edge << " as " << name << ".";
    throw (Exception(errMessage.str()));
  }
}


GlobalGraph::Edge GlobalGraph::link(Graph::NodeId nodeA, Graph::NodeId nodeB)
{
  // which ID is available?
  GlobalGraph::Edge edgeID = ++highestEdgeID_;

  // writing the new relation to the structure
  linkInNodeStructure_(nodeA, nodeB, edgeID);
  if (!directed_)
  {
    linkInNodeStructure_(nodeB, nodeA, edgeID);
  }
  linkInEdgeStructure_(nodeA, nodeB, edgeID);
  return edgeID;
}

void GlobalGraph::link(Graph::NodeId nodeA, Graph::NodeId nodeB, GlobalGraph::Edge edgeID)
{
  if (edgeStructure_.find(edgeID) != edgeStructure_.end())
    throw Exception("GlobalGraph::link : already existing edgeId " + TextTools::toString(edgeID));

  // writing the new relation to the structure
  linkInNodeStructure_(nodeA, nodeB, edgeID);
  if (!directed_)
  {
    linkInNodeStructure_(nodeB, nodeA, edgeID);
  }
  linkInEdgeStructure_(nodeA, nodeB, edgeID);
}

vector<GlobalGraph::Edge> GlobalGraph::unlink(Graph::NodeId nodeA, Graph::NodeId nodeB)
{
  // unlinking in the structure
  vector<GlobalGraph::Edge> deletedEdges; // what edges ID are affected by this unlinking
  deletedEdges.push_back(unlinkInNodeStructure_(nodeA, nodeB));

  for (vector<GlobalGraph::Edge>::iterator currEdgeToDelete = deletedEdges.begin(); currEdgeToDelete != deletedEdges.end(); currEdgeToDelete++)
  {
    unlinkInEdgeStructure_(*currEdgeToDelete);
  }

  // telling the observers
  notifyDeletedEdges(deletedEdges);

  return deletedEdges;
}

void GlobalGraph::switchNodes(Graph::NodeId nodeA, Graph::NodeId nodeB)
{
  Graph::NodeId father, son;

  nodeStructureType::iterator nodeARow = nodeStructure_.find(nodeA);
  nodeStructureType::iterator nodeBRow = nodeStructure_.find(nodeB);
  nodeStructureType::iterator nodeSonRow, nodeFatherRow;

  // Forwards
  map<GlobalGraph::Node, GlobalGraph::Edge>::iterator foundForwardRelation = nodeARow->second.first.find(nodeB);
  if (foundForwardRelation == nodeARow->second.first.end())
  {
    foundForwardRelation = nodeBRow->second.first.find(nodeA);
    if (foundForwardRelation == nodeBRow->second.first.end())
      throw Exception("GlobalGraph::exchangeNodes : no edge between nodes " + TextTools::toString(nodeA) + " and " + TextTools::toString(nodeB));
    father = nodeB;
    son = nodeA;
    nodeFatherRow = nodeBRow;
    nodeSonRow = nodeARow;
  }
  else
  {
    father = nodeA;
    son = nodeB;
    nodeFatherRow = nodeARow;
    nodeSonRow = nodeBRow;
  }

  // Edge
  GlobalGraph::Edge foundEdge = foundForwardRelation->second;

  // Backwards
  map<GlobalGraph::Node, GlobalGraph::Edge>::iterator foundBackwardsRelation = nodeSonRow->second.second.find(father);


  // Exchange
  nodeFatherRow->second.first.erase(foundForwardRelation);
  nodeSonRow->second.second.erase(foundBackwardsRelation);


  nodeSonRow->second.first[father] = foundEdge;

  nodeFatherRow->second.second[son] = foundEdge;


//    std::map<GlobalGraph::Node, std::pair<std::map<GlobalGraph::Node, GlobalGraph::Edge>, std::map<GlobalGraph::Node, GlobalGraph::Edge> > >::iterator ita = nodeStructure_.find(nodeA);

  edgeStructure_[foundEdge] = pair<Node, Node>(son, father);

  this->topologyHasChanged_();
}


GlobalGraph::Node GlobalGraph::getHighestNodeID() const
{
  return highestNodeID_;
}


GlobalGraph::Edge GlobalGraph::getHighestEdgeID() const
{
  return highestEdgeID_;
}

void GlobalGraph::unlinkInEdgeStructure_(const GlobalGraph::Edge& edge)
{
  edgeStructureType::iterator foundEdge = edgeStructure_.find(edge);
  if (foundEdge == edgeStructure_.end())
    throw Exception("GlobalGraph::unlinkInEdgeStructure_ : no edge to erase " + TextTools::toString(edge));

  edgeStructure_.erase(foundEdge);
  this->topologyHasChanged_();
}

void GlobalGraph::linkInEdgeStructure_(const GlobalGraph::Node& nodeA, const GlobalGraph::Node& nodeB, const GlobalGraph::Edge& edge)
{
  edgeStructure_[edge] = pair<Node, Node>(nodeA, nodeB);
  this->topologyHasChanged_();
}


unsigned int GlobalGraph::unlinkInNodeStructure_(const GlobalGraph::Node& nodeA, const GlobalGraph::Node& nodeB)
{
  // Forward
  nodeStructureType::iterator nodeARow = nodeStructure_.find(nodeA);
  map<GlobalGraph::Node, GlobalGraph::Edge>::iterator foundForwardRelation = nodeARow->second.first.find(nodeB);
  if (foundForwardRelation == nodeARow->second.first.end())
    throw Exception("GlobalGraph::unlinkInNodeStructure_ : no edge to erase " + TextTools::toString(nodeA) + "->" + TextTools::toString(nodeB));

  GlobalGraph::Edge foundEdge = foundForwardRelation->second;
  nodeARow->second.first.erase(foundForwardRelation);

  // Backwards
  nodeStructureType::iterator nodeBRow = nodeStructure_.find(nodeB);
  map<GlobalGraph::Node, GlobalGraph::Edge>::iterator foundBackwardsRelation = nodeBRow->second.second.find(nodeA);
  if (foundBackwardsRelation == nodeBRow->second.first.end())
    throw Exception("GlobalGraph::unlinkInNodeStructure_ : no edge to erase " + TextTools::toString(nodeB) + "<-" + TextTools::toString(nodeA));

  nodeBRow->second.second.erase(foundBackwardsRelation);

  this->topologyHasChanged_();
  return foundEdge;
}

void GlobalGraph::linkInNodeStructure_(const GlobalGraph::Node& nodeA, const GlobalGraph::Node& nodeB, const GlobalGraph::Edge& edge)
{
  std::map<GlobalGraph::Node, std::pair<std::map<GlobalGraph::Node, GlobalGraph::Edge>, std::map<GlobalGraph::Node, GlobalGraph::Edge> > >::iterator ita = nodeStructure_.find(nodeA);

  if (ita != nodeStructure_.end())
    ita->second.first.insert( pair<GlobalGraph::Node, GlobalGraph::Edge>(nodeB, edge));

  std::map<GlobalGraph::Node, std::pair<std::map<GlobalGraph::Node, GlobalGraph::Edge>, std::map<GlobalGraph::Node, GlobalGraph::Edge> > >::iterator itb = nodeStructure_.find(nodeB);

  if (itb != nodeStructure_.end())
    nodeStructure_.find(nodeB)->second.second.insert( pair<GlobalGraph::Node, GlobalGraph::Edge>(nodeA, edge));

  this->topologyHasChanged_();
}

Graph::NodeId GlobalGraph::createNode()
{
  GlobalGraph::Node newNode = highestNodeID_++;
  nodeStructure_[newNode] = std::pair<std::map<GlobalGraph::Node, GlobalGraph::Edge>, std::map<GlobalGraph::Node, GlobalGraph::Edge> >();
  this->topologyHasChanged_();

  return newNode;
}

Graph::NodeId GlobalGraph::createNodeFromNode(Graph::NodeId origin)
{
  Graph::NodeId newNode = createNode();
  link(origin, newNode);
  this->topologyHasChanged_();
  return newNode;
}

Graph::NodeId GlobalGraph::createNodeOnEdge(Graph::EdgeId edge)
{
  // origin must be an existing edge
  edgeMustExist_(edge, "");

  Graph::NodeId newNode = createNode();

  // determining the nodes on the border of the edge
  pair<GlobalGraph::Node, GlobalGraph::Node> nodes = edgeStructure_[edge];
  GlobalGraph::Node nodeA = nodes.first;
  GlobalGraph::Node nodeB = nodes.second;

  unlink(nodeA, nodeB);
  link(nodeA, newNode);
  link(newNode, nodeB);
  this->topologyHasChanged_();
  return newNode;
}


Graph::NodeId GlobalGraph::createNodeFromEdge(Graph::NodeId origin)
{
  // origin must be an existing edge
  edgeMustExist_(origin, "origin edge");

  // splitting the edge
  Graph::NodeId anchor = createNodeOnEdge(origin);

  Graph::NodeId newNode = createNodeFromNode(anchor);
  this->topologyHasChanged_();
  return newNode;
}

/*********************************************/

void GlobalGraph::registerObserver(GraphObserver* observer)
{
  if (!observers_.insert(observer).second)
    throw (Exception("This GraphObserver was already an observer of this Graph"));
  ;
}

void GlobalGraph::unregisterObserver(GraphObserver* observer)
{
  if (!observers_.erase(observer))
    throw (Exception("This GraphObserver was not an observer of this Graph"));
}


/**********************************************/

std::vector< GlobalGraph::Node > GlobalGraph::getNeighbors_(const GlobalGraph::Node& node, bool outgoing) const
{
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if (foundNode == nodeStructure_.end())
    throw (Exception("The requested node is not in the structure."));
  const std::map<GlobalGraph::Node, GlobalGraph::Edge>& forOrBack = (outgoing ? foundNode->second.first : foundNode->second.second);
  vector<GlobalGraph::Node> result;
  for (map<GlobalGraph::Node, GlobalGraph::Edge>::const_iterator currNeighbor = forOrBack.begin(); currNeighbor != forOrBack.end(); currNeighbor++)
  {
    result.push_back(currNeighbor->first);
  }

  return result;
}

std::vector< GlobalGraph::Edge > GlobalGraph::getEdges_(const GlobalGraph::Node& node, bool outgoing) const
{
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if (foundNode == nodeStructure_.end())
    throw (Exception("The requested node is not in the structure."));
  const std::map<GlobalGraph::Node, GlobalGraph::Edge>& forOrBack = (outgoing ? foundNode->second.first : foundNode->second.second);
  vector<GlobalGraph::Edge> result;
  for (map<GlobalGraph::Node, GlobalGraph::Edge>::const_iterator currNeighbor = forOrBack.begin(); currNeighbor != forOrBack.end(); currNeighbor++)
  {
    result.push_back(currNeighbor->second);
  }

  return result;
}

vector< Graph::NodeId > GlobalGraph::getIncomingNeighbors(Graph::NodeId node) const
{
  return getNeighbors_(node, false);
}

vector< Graph::EdgeId > GlobalGraph::getIncomingEdges(const Graph::NodeId node) const
{
  return getEdges_(node, false);
}

vector< Graph::NodeId > GlobalGraph::getOutgoingNeighbors(const Graph::NodeId node) const
{
  return getNeighbors_(node, true);
}

vector< Graph::EdgeId > GlobalGraph::getOutgoingEdges(const Graph::NodeId node) const
{
  return getEdges_(node, true);
}

std::unique_ptr<Graph::NodeIterator> GlobalGraph::allNodesIterator()
{
  return std::unique_ptr<Graph::NodeIterator>(new NodesIteratorClass<Graph::ALLGRAPHITER, false>(*this));
}

std::unique_ptr<Graph::NodeIterator> GlobalGraph::allNodesIterator() const
{
  return std::unique_ptr<Graph::NodeIterator>(new NodesIteratorClass<Graph::ALLGRAPHITER, true>(*this));
}


std::unique_ptr<Graph::NodeIterator> GlobalGraph::outgoingNeighborNodesIterator(Graph::NodeId node)
{
  return std::unique_ptr<Graph::NodeIterator>(new NodesIteratorClass<Graph::OUTGOINGNEIGHBORITER, false>(*this, node));
}


std::unique_ptr<Graph::NodeIterator> GlobalGraph::incomingNeighborNodesIterator(Graph::NodeId node)
{
  return std::unique_ptr<Graph::NodeIterator>(new NodesIteratorClass<Graph::INCOMINGNEIGHBORITER, false>(*this, node));
}

std::unique_ptr<Graph::NodeIterator> GlobalGraph::outgoingNeighborNodesIterator(Graph::NodeId node) const
{
  return std::unique_ptr<Graph::NodeIterator>(new NodesIteratorClass<Graph::OUTGOINGNEIGHBORITER, true>(*this, node));
}


std::unique_ptr<Graph::NodeIterator> GlobalGraph::incomingNeighborNodesIterator(Graph::NodeId node) const
{
  return std::unique_ptr<Graph::NodeIterator>(new NodesIteratorClass<Graph::INCOMINGNEIGHBORITER, true>(*this, node));
}


size_t GlobalGraph::getNumberOfNodes() const
{
  return nodeStructure_.size();
}


size_t GlobalGraph::getNumberOfEdges() const
{
  return edgeStructure_.size();
}


size_t GlobalGraph::getDegree(const Graph::NodeId node) const
{
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if (foundNode == nodeStructure_.end())
    throw Exception("GlobalGraph::getDegree : Node " + TextTools::toString(node) + " does not exist.");

  return isDirected() ? foundNode->second.first.size() + foundNode->second.second.size() : foundNode->second.first.size();
}


bool GlobalGraph::isLeaf(const Graph::NodeId node) const
{
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if (foundNode == nodeStructure_.end())
    throw Exception("GlobalGraph::isLeaf : Node " + TextTools::toString(node) + " does not exist.");

  return (!isDirected() && (foundNode->second.first.size() <= 1))
         || (isDirected() && (
               (foundNode->second.first.size() + foundNode->second.second.size() <= 1)
               || (foundNode->second.first.size() == 1 &&  foundNode->second.second.size() == 1 && foundNode->second.first.begin()->first == foundNode->second.second.begin()->first)));
}


size_t GlobalGraph::getNumberOfNeighbors(const Graph::NodeId node) const
{
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if (foundNode == nodeStructure_.end())
    throw (Exception("The requested node is not in the structure."));

  if (isDirected())
    return foundNode->second.first.size() + foundNode->second.second.size();
  else
    return foundNode->second.first.size();
}

size_t GlobalGraph::getNumberOfOutgoingNeighbors(const Graph::NodeId node) const
{
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if (foundNode == nodeStructure_.end())
    throw (Exception("The requested node is not in the structure."));
  return foundNode->second.first.size();
}

size_t GlobalGraph::getNumberOfIncomingNeighbors(const Graph::NodeId node) const
{
  nodeStructureType::const_iterator foundNode = nodeStructure_.find(node);
  if (foundNode == nodeStructure_.end())
    throw (Exception("The requested node is not in the structure."));
  return foundNode->second.second.size();
}

vector< Graph::NodeId > GlobalGraph::getNeighbors(const Graph::NodeId node) const
{
  vector<Graph::NodeId> result;
  vector<Graph::NodeId> neighborsToInsert;
  neighborsToInsert = getNeighbors_(node, false);
  result.insert(result.end(), neighborsToInsert.begin(), neighborsToInsert.end());
  neighborsToInsert = getNeighbors_(node, true);
  result.insert(result.end(), neighborsToInsert.begin(), neighborsToInsert.end());
  return result;
}

std::pair<Graph::NodeId, Graph::NodeId> GlobalGraph::getNodes(Graph::EdgeId edge) const
{
  edgeMustExist_(edge);
  edgeStructureType::const_iterator found = edgeStructure_.find(edge);
  // TODO Except if not found
  return found->second;
}

Graph::NodeId GlobalGraph::getTop(Graph::EdgeId edge) const
{
  return std::get<0>(getNodes(edge));
}

Graph::NodeId GlobalGraph::getBottom(Graph::EdgeId edge) const
{
  return std::get<1>(getNodes(edge));
}

void GlobalGraph::deleteNode(Graph::NodeId node)
{
  // checking the node
  nodeMustExist_(node, "node to delete");
  isolate_(node);

  nodeStructureType::iterator found = nodeStructure_.find(node);
  if (found == nodeStructure_.end())
    throw Exception("GlobalGraph::deleteNode : no node to erase " + TextTools::toString(node));

  nodeStructure_.erase(found);

  this->topologyHasChanged_();
}

void GlobalGraph::isolate_(GlobalGraph::Node& node)
{
  vector<Graph::NodeId> oneighbors = getOutgoingNeighbors(node);
  for (vector<Graph::NodeId>::iterator currNeighbor = oneighbors.begin(); currNeighbor != oneighbors.end(); currNeighbor++)
  {
    unlink(node, *currNeighbor);
  }

  vector<Graph::NodeId> ineighbors = getIncomingNeighbors(node);
  for (vector<Graph::NodeId>::iterator currNeighbor = ineighbors.begin(); currNeighbor != ineighbors.end(); currNeighbor++)
  {
    unlink(*currNeighbor, node);
  }
}

vector<Graph::EdgeId> GlobalGraph::getAllEdges() const
{
  vector<Graph::EdgeId> listOfEdges;
  for (edgeStructureType::const_iterator it = edgeStructure_.begin(); it != edgeStructure_.end(); it++)
  {
    listOfEdges.push_back(it->first);
  }

  return listOfEdges;
}

Graph::EdgeId GlobalGraph::getAnyEdge(Graph::NodeId nodeA, Graph::NodeId nodeB) const
{
  try
  {
    // trying in the given order A->B
    return getEdge(nodeA, nodeB);
  }
  catch (Exception e)
  {
    // didn’t work, hence trying in the opposite order B->A
    return getEdge(nodeB, nodeA);
  }
}

vector<Graph::NodeId> GlobalGraph::getAllLeaves() const
{
  vector<Graph::NodeId> listOfLeaves;
  for (nodeStructureType::const_iterator it = nodeStructure_.begin(); it != nodeStructure_.end(); it++)
  {
    if (this->isLeaf(it->first))
      listOfLeaves.push_back(it->first);
  }

  return listOfLeaves;
}

set<Graph::NodeId> GlobalGraph::getSetOfAllLeaves() const
{
  set<Graph::NodeId> listOfLeaves;
  for (nodeStructureType::const_iterator it = nodeStructure_.begin(); it != nodeStructure_.end(); it++)
  {
    if (this->isLeaf(it->first))
      listOfLeaves.insert(it->first);
  }

  return listOfLeaves;
}

vector<Graph::NodeId> GlobalGraph::getAllNodes() const
{
  vector<Graph::NodeId> listOfNodes;
  for (nodeStructureType::const_iterator it = nodeStructure_.begin(); it != nodeStructure_.end(); it++)
  {
    listOfNodes.push_back(it->first);
  }

  return listOfNodes;
}

vector<Graph::NodeId> GlobalGraph::getAllInnerNodes() const
{
  vector<Graph::NodeId> listOfInNodes;
  for (nodeStructureType::const_iterator it = nodeStructure_.begin(); it != nodeStructure_.end(); it++)
  {
    if (this->getDegree(it->first) >= 2)
      listOfInNodes.push_back(it->first);
  }

  return listOfInNodes;
}


void GlobalGraph::fillListOfLeaves_(const GlobalGraph::Node& startingNode, vector<GlobalGraph::Node>& foundLeaves, const GlobalGraph::Node& originNode, unsigned int maxRecursions) const
{
  const vector<Graph::NodeId> neighbors = getNeighbors(startingNode);
  if (neighbors.size() > 1)
  {
    if (maxRecursions > 0)
      for (vector<Node>::const_iterator currNeighbor = neighbors.begin(); currNeighbor != neighbors.end(); currNeighbor++)
      {
        if (*currNeighbor != originNode)
          fillListOfLeaves_(*currNeighbor, foundLeaves, startingNode, maxRecursions - 1);
      }
  }
  else
  {
    foundLeaves.push_back(startingNode);
  }
}


std::vector<Graph::NodeId> GlobalGraph::getLeavesFromNode(const Graph::NodeId node, unsigned int maxDepth) const
{
  vector<Graph::NodeId> listOfLeaves;
  fillListOfLeaves_(node, listOfLeaves, node, maxDepth);
  return listOfLeaves;
}

void GlobalGraph::nodeToDot_(const GlobalGraph::Node& node, ostream& out,  std::set<std::pair<Node, Node> >& alreadyFigured) const
{
  bool theEnd = true;
  const std::map<Node, Edge>& children = nodeStructure_.at(node).first;
  for (map<Node, Edge>::const_iterator currChild = children.begin(); currChild != children.end(); currChild++)
  {
    if (alreadyFigured.find(pair<Node, Node>(node, currChild->first)) != alreadyFigured.end() || (!directed_ && alreadyFigured.find(pair<Node, Node>(currChild->first, node)) != alreadyFigured.end()))
      continue;
    alreadyFigured.insert(pair<Node, Node>(node, currChild->first));
    theEnd = false;
    out << node << (directed_ ? " -> " : " -- ");
    nodeToDot_(currChild->first, out, alreadyFigured);
  }

  const std::map<Node, Edge>& fathers = nodeStructure_.at(node).second;
  for (map<Node, Edge>::const_iterator currFath = fathers.begin(); currFath != fathers.end(); currFath++)
  {
    if (alreadyFigured.find(pair<Node, Node>(currFath->first, node)) != alreadyFigured.end() || (!directed_ && alreadyFigured.find(pair<Node, Node>(node, currFath->first)) != alreadyFigured.end()))
      continue;
    alreadyFigured.insert(pair<Node, Node>(currFath->first, node));
    theEnd = false;
    out << node << (directed_ ? " <- " : " -- ");
    nodeToDot_(currFath->first, out, alreadyFigured);
  }
  if (theEnd)
    out << node << ";\n    ";
}

bool GlobalGraph::isTree() const
{
  set<GlobalGraph::Node> metNodes;
  bool nodesAreMetOnlyOnce = nodesAreMetOnlyOnce_(root_, metNodes, root_);

  if (!nodesAreMetOnlyOnce)
    return false;
  // now they have only been met at most once, they have to be met at least once
  bool noNodeMissing = true;
  for (nodeStructureType::const_iterator currNode = nodeStructure_.begin(); noNodeMissing && currNode != nodeStructure_.end(); currNode++)
  {
    noNodeMissing = (metNodes.find(currNode->first) != metNodes.end());
  }
  return noNodeMissing;
}


bool GlobalGraph::nodesAreMetOnlyOnce_(const GlobalGraph::Node& node, set< GlobalGraph::Node >& metNodes, const GlobalGraph::Node& originNode) const
{
  // insert().second <=> not yet in the set
  bool neverMetANodeMoreThanOnce = metNodes.insert(node).second;
  vector<Graph::NodeId> neighbors = getOutgoingNeighbors(node);
  for (vector<Graph::NodeId>::iterator currNeighbor = neighbors.begin(); neverMetANodeMoreThanOnce && currNeighbor != neighbors.end(); currNeighbor++)
  {
    if (*currNeighbor == originNode)
      continue;
    neverMetANodeMoreThanOnce = nodesAreMetOnlyOnce_(*currNeighbor, metNodes, node);
  }
  return neverMetANodeMoreThanOnce;
}

bool GlobalGraph::isDA() const
{
  GlobalGraph gg(*this);

  // Algo: remove recursively all nodes with no sons from graph

  std::vector<Graph::NodeId> vL;

  std::unique_ptr<Graph::NodeIterator> it = gg.allNodesIterator();
  for ( ; !it->end(); it->next())
  {
    if (gg.getNumberOfOutgoingNeighbors(**it) == 0)
      vL.push_back(**it);
  }

  while (vL.size() != 0)
  {
    for (std::vector<Graph::NodeId>::iterator it2(vL.begin()); it2 != vL.end(); it2++)
    {
      gg.deleteNode(*it2);
    }

    if (gg.getNumberOfNodes() == 0)
      return true;

    vL.clear();

    it = gg.allNodesIterator();
    for ( ; !it->end(); it->next())
    {
      if (gg.getNumberOfOutgoingNeighbors(**it) == 0)
        vL.push_back(**it);
    }
  }

  return false;
}


void GlobalGraph::orientate()
{
  if (!isDirected())
    makeDirected();

  GlobalGraph gg(*this);

  // Algo: remove recursively all nodes from graph, starting with
  // root_

  Graph::NodeId node = root_;
  std::set<Graph::NodeId> nextNodes;
  nextNodes.insert(node);

  while (gg.getNumberOfNodes() != 0)
  {
    // look for the next node to be treated
    Graph::NodeId nbgg = 0;

    // first node with one neighbor (ie no choice on orientation)

    std::set<Graph::NodeId>::iterator it = nextNodes.begin();
    for ( ; it != nextNodes.end(); it++)
    {
      if (gg.getNumberOfNeighbors(*it) <= 1)
        break;
    }

    // if none, look for node wih minimum number of fathers
    if (it == nextNodes.end())
    {
      size_t nbF = numeric_limits<size_t>::infinity();
      it = nextNodes.begin();

      for ( ; it != nextNodes.end(); it++)
      {
        size_t nbFi = gg.getNumberOfIncomingNeighbors(*it);
        if (nbF == 0)
        {
          nbgg = *it;
          break;
        }
        else
        {
          if (nbFi < nbF)
          {
            nbgg = *it;
            nbF = nbFi;
          }
        }
      }
    }
    else
      nbgg = *it;

    // next orient edges from this node and catch neighbors
    std::vector<Graph::NodeId> vL = gg.getIncomingNeighbors(nbgg);
    for (std::vector<Graph::NodeId>::iterator it2(vL.begin()); it2 != vL.end(); it2++)
    {
      switchNodes(nbgg, *it2);

      nextNodes.insert(*it2);
    }

    vL = gg.getOutgoingNeighbors(nbgg);
    for (std::vector<Graph::NodeId>::iterator it2(vL.begin()); it2 != vL.end(); it2++)
    {
      nextNodes.insert(*it2);
    }

    gg.deleteNode(nbgg);
    nextNodes.erase(nbgg);
  }
}

void GlobalGraph::setRoot(Graph::NodeId newRoot)
{
  nodeMustExist_(newRoot, "new root");
  root_ = newRoot;
}

Graph::NodeId GlobalGraph::getRoot() const
{
  return root_;
}


bool GlobalGraph::isDirected() const
{
  return directed_;
}

void GlobalGraph::makeDirected()
{
  if (directed_)
    return;
  // save and clean the undirectedStructure
  nodeStructureType undirectedStructure = nodeStructure_;
  for (nodeStructureType::iterator it = nodeStructure_.begin(); it != nodeStructure_.end(); it++)
  {
    it->second = std::pair<std::map<Node, Edge>, std::map<Node, Edge> >();
  }
  // copy each relation once, without the reciprocal link
  // (first met, first kept)
  // eg: A - B in undirected is represented as A->B and B->A
  //     in directed, becomes A->B only
  std::set<pair<Node, Node> > alreadyConvertedRelations;
  for (nodeStructureType::iterator currNodeRow = undirectedStructure.begin(); currNodeRow != undirectedStructure.end(); currNodeRow++)
  {
    Node nodeA = currNodeRow->first;

    for (map<Node, Edge>::iterator currRelation = currNodeRow->second.first.begin(); currRelation != currNodeRow->second.first.end(); currRelation++)
    {
      Node nodeB = currRelation->first;
      Edge edge = currRelation->second;
      if (alreadyConvertedRelations.insert(pair<Node, Node>(min(nodeA, nodeB), max(nodeA, nodeB))).second)
        linkInNodeStructure_(nodeA, nodeB, edge);
    }
  }
  directed_ = true;
  this->topologyHasChanged_();
}

void GlobalGraph::makeUndirected()
{
  if (!directed_)
    return;
  if (containsReciprocalRelations())
    throw Exception("Cannot make an undirected graph from a directed one containing reciprocal relations.");
  // save and clean the undirectedStructure
  nodeStructureType directedStructure = nodeStructure_;
  for (nodeStructureType::iterator it = nodeStructure_.begin(); it != nodeStructure_.end(); it++)
  {
    it->second = std::pair<std::map<Node, Edge>, std::map<Node, Edge> >();
  }
  // copy each relation twice, making the reciprocal link
  // eg: A - B in directed is represented as A->B
  //     in undirected, becomes A->B and B->A
  for (nodeStructureType::iterator currNodeRow = directedStructure.begin(); currNodeRow != directedStructure.end(); currNodeRow++)
  {
    Node nodeA = currNodeRow->first;
    for (map<Node, Edge>::iterator currRelation = currNodeRow->second.first.begin(); currRelation != currNodeRow->second.first.end(); currRelation++)
    {
      Node nodeB = currRelation->first;
      Edge edge = currRelation->second;
      linkInNodeStructure_(nodeA, nodeB, edge);
      linkInNodeStructure_(nodeB, nodeA, edge);
    }
  }
  directed_ = false;
  this->topologyHasChanged_();
}

bool GlobalGraph::containsReciprocalRelations() const
{
  if (!directed_)
    throw Exception("Cannot state reciprocal link in an undirected graph.");
  std::set<pair<Node, Node> > alreadyMetRelations;
  for (nodeStructureType::const_iterator currNodeRow = nodeStructure_.begin(); currNodeRow != nodeStructure_.end(); currNodeRow++)
  {
    Node nodeA = currNodeRow->first;
    for (map<Node, Edge>::const_iterator currRelation = currNodeRow->second.first.begin(); currRelation != currNodeRow->second.first.end(); currRelation++)
    {
      Node nodeB = currRelation->first;
      if (!alreadyMetRelations.insert(pair<Node, Node>(min(nodeA, nodeB), max(nodeA, nodeB))).second)
        return true;
    }
  }
  return false;
}

std::unique_ptr<Graph::EdgeIterator> GlobalGraph::allEdgesIterator()
{
  return std::unique_ptr<Graph::EdgeIterator>(new EdgesIteratorClass<Graph::ALLGRAPHITER, false>(*this));
}

std::unique_ptr<Graph::EdgeIterator> GlobalGraph::outgoingEdgesIterator(Graph::NodeId node)
{
  return std::unique_ptr<Graph::EdgeIterator>(new EdgesIteratorClass<Graph::OUTGOINGNEIGHBORITER, false>(*this, node));
}

std::unique_ptr<Graph::EdgeIterator> GlobalGraph::incomingEdgesIterator(Graph::NodeId node)
{
  return std::unique_ptr<Graph::EdgeIterator>(new EdgesIteratorClass<Graph::INCOMINGNEIGHBORITER, false>(*this, node));
}

std::unique_ptr<Graph::EdgeIterator> GlobalGraph::allEdgesIterator() const
{
  return std::unique_ptr<Graph::EdgeIterator>(new EdgesIteratorClass<Graph::ALLGRAPHITER, true>(*this));
}

std::unique_ptr<Graph::EdgeIterator> GlobalGraph::outgoingEdgesIterator(Graph::NodeId node) const
{
  return std::unique_ptr<Graph::EdgeIterator>(new EdgesIteratorClass<Graph::OUTGOINGNEIGHBORITER, true>(*this, node));
}

std::unique_ptr<Graph::EdgeIterator> GlobalGraph::incomingEdgesIterator(Graph::NodeId node) const
{
  return std::unique_ptr<Graph::EdgeIterator>(new EdgesIteratorClass<Graph::INCOMINGNEIGHBORITER, true>(*this, node));
}

Graph::EdgeId GlobalGraph::getEdge(Graph::NodeId nodeA, Graph::NodeId nodeB) const
{
  nodeStructureType::const_iterator firstNodeFound = nodeStructure_.find(nodeA);
  if (firstNodeFound == nodeStructure_.end())
    throw (Exception("The fist node was not the origin of an edge."));
  map<Node, Edge>::const_iterator secondNodeFound = firstNodeFound->second.first.find(nodeB);
  if (secondNodeFound == firstNodeFound->second.first.end())
    throw (Exception("The second node was not in a relation with the first one."));
  return secondNodeFound->second;
}

vector<Graph::EdgeId> GlobalGraph::getEdges(Graph::NodeId node) const
{
  vector<Graph::EdgeId> result;
  vector<Graph::EdgeId> edgesToInsert;
  edgesToInsert = getEdges_(node, false);
  result.insert(result.end(), edgesToInsert.begin(), edgesToInsert.end());
  edgesToInsert = getEdges_(node, true);
  result.insert(result.end(), edgesToInsert.begin(), edgesToInsert.end());
  return result;
}

void GlobalGraph::outputToDot(ostream& out, const std::string& name) const
{
  out << (directed_ ? "digraph" : "graph") << " " << name << " {\n    ";
  set<pair<Node, Node> > alreadyFigured;
  nodeToDot_(root_, out, alreadyFigured);
  out << "\r}" << endl;
}

void GlobalGraph::notifyDeletedEdges(const vector<Graph::EdgeId>& edgesToDelete) const
{
  for (set<GraphObserver*>::iterator currObserver = observers_.begin(); currObserver != observers_.end(); currObserver++)
  {
    (*currObserver)->deletedEdgesUpdate(edgesToDelete);
  }
}

void GlobalGraph::notifyDeletedNodes(const vector<Graph::NodeId>& nodesToDelete) const
{
  for (set<GraphObserver*>::iterator currObserver = observers_.begin(); currObserver != observers_.end(); currObserver++)
  {
    (*currObserver)->deletedNodesUpdate(nodesToDelete);
  }
}
