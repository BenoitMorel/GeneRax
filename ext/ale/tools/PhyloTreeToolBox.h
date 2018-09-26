// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 12/01/17.
//

#ifndef PHYLASOLVER_PHYLOTREETOOLBOX_H
#define PHYLASOLVER_PHYLOTREETOOLBOX_H

//Includes Bpp
#include <Bpp/Phyl/Tree/PhyloTree.h>

//Includes standards
#include <string>
#include <unordered_set>
#include <stack>
#include <queue>
#include <cmath>

#include <ale/Defines.h>

//Includes containers
#include <ale/containers/Cost.h>
#include <ale/containers/GeneMap.h>

class PhyloTreeToolBox {
  /*!
   * \class PhyloTreeToolBox
   * \brief PhyloTreeToolBox provides methods to explore and modify a bpp::PhyloTree as pruning, Post-order traversal edition, bpp::PhyloTree cloning, etc.
   * \details
   *   Get all nodes according to a post-order traversal:
   *
   *        auto nodes = PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(my_PhyloTree);
   *
   *   Get all nodes according to a name pattern:
   *
   *        auto nodes_selected = PhyloTreeToolBox::getNodesFromNamePattern(nodes, "homo_*");
   *
   */

private:
  /************************************
   * Private methods
   */


  /// \brief Recursive function for in-order traversal.
  /// \return Nothing, used to fill a list of nodes according to a visiting Node node and a Condition condition.
  template<typename Condition>
  static void _fillInOrderRecursive_(const bpp::PhyloTree &tree /// tree to explore.
      , const Node &node /// node which is visited.
      , std::list<Node> &nodes /// list of nodes in in-order traversal.
      , const Condition &condition /// condition to push back the node in the resulting list.
  ) {
    bool isLeaf = tree.isLeaf(node);

    std::vector<Node> sons;

    if(not isLeaf) {
      sons = tree.getSons(node);
      _fillInOrderRecursive_(tree, sons.front(), nodes, condition);
    }

    if(condition(node))
      nodes.push_back(node);

    if(sons.size() > 1) {
      for(std::size_t i = 1; i < sons.size(); ++i)
        _fillInOrderRecursive_(tree, sons.at(i), nodes, condition);
    }
  }

  /// \brief Recursive function for pre-order traversal.
  /// \return Nothing, used to fill a list of nodes according to a visiting Node node and a Condition condition.
  template<typename Condition>
  static void _fillInPreOrderRecursive_(const bpp::PhyloTree &tree /// tree to explore.
      , const Node &node /// node which is visited.
      , std::list<Node> &nodes /// list of nodes in pre-order traversal.
      , const Condition &condition /// condition to push back the node in the resulting list.
  ) {
    if(condition(node))
      nodes.push_back(node);

    if(not tree.isLeaf(node))
      for(auto& son: tree.getSons(node))
        _fillInPreOrderRecursive_(tree, son, nodes, condition);
  }

  /// \brief Recursive function for post-order traversal.
  /// \return Nothing, used to fill a list of nodes according to a visiting Node node and a Condition condition.
  template<typename Condition>
  static void _fillInPostOrderRecursive_(const bpp::PhyloTree &tree /// tree to explore.
      , const Node &node /// node which is visited.
      , std::list<Node> &nodes /// list of nodes in post-order traversal.
      , const Condition& condition /// condition to push back the node in the resulting list.
  ) {
    if(not tree.isLeaf(node))
      for(auto& son: tree.getSons(node))
        _fillInPostOrderRecursive_(tree, son, nodes, condition);

    if(condition(node))
      nodes.push_back(node);
  }


  /// \brief Recursive tree cloning. Clone only tree topology: clone and original tree shared the same nodes.
  static void _recursiveTreeCloning_(
      const bpp::PhyloTree &tree /// Tree to clone.
      , const std::shared_ptr<bpp::PhyloNode> &current_node
      , const std::shared_ptr<bpp::PhyloNode> &father_node
      , bpp::PhyloTree &newtree
      , const bool keep_branch_lengths
      , std::unordered_set<std::shared_ptr<bpp::PhyloNode>>& created
      , const bool verbose
  );

  /// \brief Returns a list of node (in Post-order traversal) according to a "push_back" condition and an other to stop exploration. Iterative version.
  /// \return A list of nodes according to a post-order traversal, a condition and a subtree root.
  template<typename Get_condition, typename Stop_condition>
  static std::list<Node> PostOrderTraversalIterative_list(
      const bpp::PhyloTree &tree /// PhyloTree to visit.
      , const Get_condition &get_condition /// Getter condition. Takes a Node in parameter and returns a boolean. If the condition returns true, the node is retained.
      , const Stop_condition &stop_condition /// Stop condition. Takes a Node in parameter and returns a boolean. If the condition return true, the traversal is stopped.
      , const Node &root = nullptr /// Node to start the traversal, default is the root of the tree.
  ) {
    std::stack<Node> s1;
    std::list<Node> s2; //used as a second stack which contains nodes in the right order.
    Node current = (root) ? root : tree.getRoot();
    s1.push(current);
    while (not s1.empty()) {
      current = s1.top();
      s1.pop();
      if (get_condition(current))
        s2.push_front(current);

      if (stop_condition(current))
        return s2;

      if (not tree.isLeaf(current))
        for (auto son: tree.getSons(current))
          s1.push(son);
    }

    return s2;
  }

protected:

public:
  //=========================================================
  // SETTERS, transform a tree, modify, remove and add nodes.
  //=========================================================

  /// \brief Apply a function (from lambda or functor) to nodes in a given tree according to a post-order traversal.
  /// \return Nothing but modify nodes in tree.
  template<typename Func>
  static void applyInPOT(bpp::PhyloTree &tree, const Func &fun, const Node &node = nullptr) {
    Node root = node ? node : tree.getRoot();
    auto nodes = getNodesInPostOrderTraversalIterative_list(tree, root);
    std::for_each(nodes.begin(), nodes.end(), fun);
  }

  /// \brief Add father node and its sons in a tree. Be careful nodes needs to be created before.
  static void addNodeWithSonsInPostOrder(bpp::PhyloTree &tree /// Tree to modify.
      , const std::shared_ptr<bpp::PhyloNode> &father /// Father node.
      , const std::shared_ptr<bpp::PhyloNode> &son_left /// Son left of the father.
      , const std::shared_ptr<bpp::PhyloNode> &son_right /// Son right of the father.
      , const std::shared_ptr<bpp::PhyloBranch> &son_left_branch = 0 /// Branch between father and son_left.
      , const std::shared_ptr<bpp::PhyloBranch> &son_right_branch = 0 /// Branch between fahter and son_right.
      , const bool verbose = false
  );

  /// \brief Contract tree branches according to a maximal threshold of branch support. Branch with a support below the
  ///        threshold are contracted.
  /// \return Returns a list of nodes which are removed from the tree.
  static std::list<Node> mergeWeakBranches(bpp::PhyloTree &tree ///tree to contract
      , const double threshold ///minimum value of branch support to keep.
      , const bool inferior_only = DEFAULT_THRESHOLD_COMPARISON_INFERIOR_ONLY ///do not contract branches if branch support equal to the threshold.
  );

  /// \brief Prune tree of artificial genes.
  static void removeArtificialGenes(bpp::PhyloTree &tree);

  /// \brief Removes nodes of a tree matching with a specific given condition.
  template<typename Func>
  static void removeNodesInPostOrderTraversal(
      bpp::PhyloTree &tree /// Tree to modify.
      , const Node &root /// Root of the subtree.
      , const Func &condition /// Condition to remove a node.
      , const bool verbose = false
  ) {
    PhyloTreeToolBox::applyInPOT(tree, [&tree, &condition, &verbose](const Node& node){
      if(condition(node, tree))
        PhyloTreeToolBox::removeNode(tree, node, verbose);
    });
  }

  /// \brief Removes a node from a given tree.
  /// \return Nothing but modify the tree by removing a node.
  static void removeNode(
      bpp::PhyloTree &tree /// Tree to modify.
      , const Node &node /// Node to remove.
      , const bool verbose = false
  );

  /// \brief Reset Node and edges ids in a bpp::Phylotree.
  /// \return Nothing but changes Node Ids in tree.
  static void resetNodeIdInPostOrder(bpp::PhyloTree &tree);

  /// \brief Remove branch length information.
  /// \return Nothing but changes tree.
  static void removeBranchLengths(bpp::PhyloTree &tree) {
    auto edges = tree.getAllEdges();
    for(auto& edge: edges) {
      if(edge->hasLength()) {
        edge->deleteLength();
      }
    }
  };


  /* Simplifying trees, pruning*/
  /// \brief Prune a species tree of leaves which are not specified in the SpeciesGeneMap.
  /// \return Nothing but edit the tree.
  static void pruneTree_v0(bpp::PhyloTree &speciestree,
                    const SpeciesGeneMap &genemap,
                    const bool verbose = false);

  //=========================================================
  // GETTERS
  // Get or find nodes according to conditions or order, get node properties, generate distance matrix, clone tree, etc.
  //=========================================================

  /// \brief Get Nodes (according to a particular given condition) in order, returned in a list. Recursive version.
  /// \return List of nodes according to an in-order traversal.
  template<typename Condition>
  static std::list<Node> InOrderTraversalRecursive_list(
      const bpp::PhyloTree& tree /// Tree to retrieve nodes.
      , const Condition& condition /// Condition to append node in the resulting list.
      , const Node& root = nullptr /// Root of the   subtree, default to the root of the entire tree.
  ) {
    std::list<Node> nodes;
    _fillInOrderRecursive_(tree, (root ? root : tree.getRoot()), nodes, condition);
    return nodes;
  }

  /// \brief Get Nodes (according to a particular given condition) in pre order, returned in a list. Recursive version.
  /// \return List of nodes according to a pre-order traversal.
  template<typename Condition>
  static std::list<Node>
  PreOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree /// Tree to retrieve nodes.
      , const Condition& condition /// Condition to append node in the resulting list.
      , const Node &root = nullptr /// Root of the subtree, default to the root of the entire tree.
  ) {
    std::list<Node> nodes;
    _fillInPreOrderRecursive_(tree, root ? root : tree.getRoot(), nodes, condition);
    return nodes;
  }

  /// \brief Get Nodes (according to a particular given condition) in post order, returned in a list. Recursive version.
  /// \return List of nodes according to a post-order traversal.
  template<typename Condition>
  static std::list<Node>
  PostOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree /// Tree to retrieve nodes.
      , const Condition& condition /// Condition to append node in the resulting list.
      , const Node &root = nullptr /// Root of the subtree, default to the root of the entire tree.
  ) {
    std::list<Node> nodes;
    _fillInPostOrderRecursive_(tree, root ? root : tree.getRoot(), nodes, condition);
    return nodes;
  }

  /// \brief Get nodes in order traversal. Recursive version.
  /// \return List of all nodes of a subtree according to an in-order traversal.
  static std::list<Node>
  getNodesInOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree /// Tree to explore.
      , const Node &root = nullptr /// Root of the subtree, default to the root of the entire tree.
  );

  /// \brief Get nodes in a pre-order traversal. Recursive version.
  /// \return List of all nodes of a subtree according to a pre-order traversal.
  static std::list<Node>
  getNodesInPreOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree /// Tree to explore.
      , const Node &root = nullptr /// Root of the subtree, default to the root of the entire tree.
  );

  /// \brief Get nodes in a post-order traversal. Recursive version.
  /// \return List of all nodes of a subtree according to a post-order traversal.
  static std::list<Node>
  getNodesInPostOrderTraversalRecursive_list(
      const bpp::PhyloTree &tree /// Tree to explore.
      , const Node &root = nullptr /// Root of the subtree, default to the root of the entire tree.
  );

  /// \brief Check if almost one node of the tree returns true according to a particular condition.
  template<typename Condition>
  static bool hasNode(
      const bpp::PhyloTree &tree /// Tree to check.
      , const Condition &condition /// Lambda condition to check, like [](const std::shared_ptr<bpp::PhyloNode>& node){...}.
      , const Node &root = nullptr /// Root of the subtree to explore, default is the root of the entire tree.
  ){
    return PostOrderTraversalIterative_list(tree, condition, condition, root).size() > 0;
  }

  /// \brief Get nodes in a post-order traversal. Iterative version.
  /// \return List of all nodes of a subtree according to a post-order traversal.
  static std::list<Node>
  getNodesInPostOrderTraversalIterative_list(
      const bpp::PhyloTree &tree /// Tree to explore.
      , const Node &root = nullptr /// Root of the subtree, default is the root of the given tree.
  );

  /// \brief Get nodes in a post-order traversal according to a specific condition.
  /// \return List of all nodes of a subtree according to a post-order traversal and which in respect of a given condition.
  template<typename Condition>
  static std::list<Node>
  getNodesInPostOrderTraversalIterative_list(
      const bpp::PhyloTree &tree /// Tree to explore.
      , const Condition &condition /// Condition to append node in the resulting list.
      , const Node &root = nullptr /// Root of the subtree to explore, default is the root of the entire tree.
  ) {
    return PostOrderTraversalIterative_list(tree, condition, [](const Node& node) { return false; }, root);
  }

  /// \brief Get names from a container of nodes using two iterators.
  /// \return Vector of names of all nodes in a container.
  template<typename Container_iterator>
  static std::vector<std::string> getNodeNames(Container_iterator begin, const Container_iterator& end){
    std::vector<std::string> node_names((unsigned long) std::distance(begin, end));

    std::generate(node_names.begin(), node_names.end(), [&begin](){
      std::string name = "";

      if((*begin)->hasName())
        name = (*begin)->getName();

      begin++;
      return name;
    });

    return node_names;
  }

  /// \brief Returns all internal nodes under a given node (default = tree's root).
  /// \return List of all internal nodes under a node.
  static std::list<Node> getInternalNodes(const bpp::PhyloTree& tree, const Node& root = nullptr);

  /// \brief Returns all internal node names under a give node (default = tree's root).
  /// \return Vector of names.
  static std::vector<std::string> getInternalNodeNames(const bpp::PhyloTree& tree, const Node& root = nullptr);

  /// \brief Returns all leaves under a node (default = tree's root).
  /// \return List of std::shared_ptr<bpp::PhyloNode>.
  static std::list<Node> getLeaves(const bpp::PhyloTree& tree, const Node& root = nullptr);

  /// \brief Get all leaves names under a node (default = tree's root).
  /// \return Vector of std::string.
  static std::vector<std::string> getLeavesNames(const bpp::PhyloTree& tree, const Node& root = nullptr);

  /// \brief Get the distance between two nodes.
  static double getDistanceBetweenTwoNodes(const Node &nodeA, const Node &nodeB, const bpp::PhyloTree &tree);

  /// \brief Check if there is a distance between two nodes.
  static bool hasDistanceBetweenTwoNodes(const Node &nodeA, const Node &nodeB, const bpp::PhyloTree& tree);

  /// \brief Sort Nodes vector by their index in the bpp::PhyloTree. Useful when nodes are associated with an id according to a
  /// post-order traversal.
  static void sortNodesByIndex(std::vector<std::shared_ptr<bpp::PhyloNode>> &nodes, const bpp::PhyloTree &tree,
                               const bool ascending = true);

  /// \brief Sort nodes according to their name.
  static void sortNodesByNameLength(std::vector<std::shared_ptr<bpp::PhyloNode>> &nodes, const bool ascending = true);

  /// \brief Returns the node which has the same name (first occurence).
  static std::shared_ptr<bpp::PhyloNode>
  getNodeFromName(const std::vector<Node> nodes, const std::string& name, const bool case_sensitive = false);

  /// \brief Returns the node which has the same name (first occurence).
  static std::shared_ptr<bpp::PhyloNode>
  getNodeFromName(const bpp::PhyloTree &tree, const std::string& name, const bool case_sensitive = false, const bool leaves_only = true);

  /// \brief Check if a node name has a match with a pattern (supports only "*" pattern).
  static bool node_name_match_with_pattern(const std::string& node_name, const std::string& pattern, const bool case_sensitive = false);

  /// \brief Returns indexes of nodes (in a std::list) with a name following the given name pattern.
  static std::list<std::size_t>
  getNodeIndexesFromNamePattern(const std::vector<Node> &nodes, const std::string& pattern, const bool case_sensitive = false);

  /// \brief Returns indexes of nodes with a name following the given name pattern.
  static std::vector<std::shared_ptr<bpp::PhyloNode>>
  getNodesFromNamePattern(const std::vector<Node> &nodes, const std::string& pattern, const bool case_sensitive = false);

  /// \brief Returns the last common ancestor (lca) of two nodes given nodes.
  static Node getCommonAncestor(const bpp::PhyloTree &tree, const Node &nodeA, const Node &nodeB);

  /// \brief Returns the last common ancestor (lca) of a set of nodes.
  static Node getLastCommonAncestor(std::vector<Node> nodes, const bpp::PhyloTree &tree);

  /// \brief Get brothers of a given node.
  static std::vector<Node> getNodeBrothers(const bpp::PhyloTree &tree, const Node &node);

  /// \brief Clone a tree topology without cloning nodes. Clone and original tree share same bpp::PhyloNode.
  static std::shared_ptr<bpp::PhyloTree> cloneTree(
      const bpp::PhyloTree &tree /// Original tree to clone.
      , const Node &root = nullptr /// Clone a subtree, according to a specific root (default is the root of the tree).
      , const bool keep_branch_lengths = true /// Keep branch lengths in resulting tree.
      , const bool verbose = false /// Verbose.
  );

  /// \brief Get branch support frequencies for a set of trees. The returned map contains in key: the branch support and the
  /// frequency as value.

  /// \brief Get branch support frequencies. The returned map contains in key: the branch support and the frequency as value.

  /// \brief Indicates if the given tree is binary (true) or not (false).
  static bool isBinary(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Compare two tree topologies. Returns true if trees are similar, and false if there is a difference.
  ///        Branch lengths are not compared.
  static bool compare_trees_with_days_algorithm(const bpp::PhyloTree &treeA, const bpp::PhyloTree &treeB,
                                         const bool verbose = false);

  /// \brief Returns true if the tree contains at least one artificial gene.
  static bool hasArtificialGenes(const bpp::PhyloTree &tree, const Node& root = nullptr);

  /// \brief Check if a node is an artificial gene. This property is defined in bpp::PhyloNode properties.
  static bool isArtificalGene(const Node &node);

  /// \brief Get Node occurrence. Occurrence is a particular integer used in Polytomysolver algorithm.
  static std::size_t getNodeOccurence(const Node &node);

  /// \brief Returns a list of nodes matching with a specified outdegree (number of sons).
  static std::list<Node> find_nodes_according_to_its_outdegree(const bpp::PhyloTree &tree, const size_t &outdegree,
                                                        const Node &root = nullptr);

  /// \brief Returns a list of multifurcation roots. Note: polytomy = multifurcation.
  /// \note We consider a polytomy (= multifurcation) as a node which has 3 or more sons.
  static std::list<Node> findPolytomies(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Returns a list of bifurcation roots.
  /// \note We consider a bifurcation as a node which has exactly two sons.
  static std::list<Node> findBifurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Returns a list of monofurcation roots. A monofurcation is an internal node with only one son.
  /// \note We consider a monofurcation as a node which has only one son.
  static std::list<Node> findMonofurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Check if a tree contains at least one monofurcation. A monofurcation is an internal node with only one son.
  static bool hasMonofurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Sort internal nodes according to their outdegree. The first element of the pair is a list of bifurcations and the
  ///       second a list of multifurcations (= polytomy).
  /// \note We consider bifurcation as a node which contains exactly two sons and a multifurcation (= polytomy) as a
  ///       node which has more than two sons.
  static std::pair<std::list<Node>, std::list<Node>>
  getFurcations(const bpp::PhyloTree &tree, const Node &root = nullptr);

  /// \brief Deduce events in a given bifurcation according to the species tree, and a map.
  /// \return A map of Events (duplication, loss) with their number.
  static std::map<Event, std::size_t> getGeneEventsInBifurcation(const Node &node /// Gene node of interest
      , const bpp::PhyloTree &genetree /// Gene tree
      , const bpp::PhyloTree &speciestree /// Species tree
      , const SpeciesGeneMap &map /// Map
  );

  /// Check if each node of a tree is connected to a father.
  static bool allEdgesAreCorrect(const bpp::PhyloTree& tree);

  /// Deduce event in a gene tree node. GeneTree needs to be reconciled.
  static Event getEvent(const bpp::PhyloTree& genetree, const SpeciesGeneMap& map, const Node& node);

};

#endif //PHYLASOLVER_PHYLOTREETOOLBOX_H
