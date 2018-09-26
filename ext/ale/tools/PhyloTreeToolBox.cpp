// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 12/01/17.
//

#include "PhyloTreeToolBox.h"

// Include standards
#include <cassert>

// Include Treerecs containers
#include <ale/containers/NodeProperty.h>

// Include Treerecs tools
#include "Utils.h"

void PhyloTreeToolBox::pruneTree_v0(bpp::PhyloTree &speciestree,
                                    const SpeciesGeneMap &genemap,
                                    const bool verbose) {
  ///
  ///\brief Cut and edit species tree branches according to the presence or not of genes in leaves.
  ///\note The SpeciesGeneMap changes the tree in parameter.
  ///

  // FIRST STEP :: prune from leaves
  if(verbose) std::cout << "Pruning tree: " << speciestree << " with " << __FUNCTION__ << std::endl;
  auto species_leafs = speciestree.getAllLeaves(); //the vector of all leaves in the species tree
  std::size_t noname_number = 0; //We can set names to anonymous nodes, so we are going to name these nodes
  while(not species_leafs.empty()){ // while the leaves to evaluate list is not empty
    std::shared_ptr<bpp::PhyloNode>& current_species = species_leafs[0]; //take the first leaf.
    if(verbose) std::cout << "\tWorking on leaf " << current_species << "..." << std::endl;
    if(not genemap.hasGenes(current_species)) { //if there is no gene associated to the species "current_species"
      if(verbose) std::cout << "\t\tLeaf '" << current_species << "' has no gene" << std::endl;
      //And, if the the sister of the leaf is empty too
      if(speciestree.hasFather(current_species)) {//if the node has a father, looking for in brothers if gene exists
        auto father = speciestree.getFather(current_species); //first get the father node
        if(father != speciestree.getRoot()) { //if current_node is not a leaf of the tree root
          std::vector<decltype(father)> sons = speciestree.getSons(father); //then get sons of father
          bool father_has_no_son_with_gene = true;
          for (auto &son: sons) {
            if (genemap.hasGenes(son) or not speciestree.isLeaf(
                son)) { //if the son is an internal node or the leaf has gene: do not delete the father -> false.
              father_has_no_son_with_gene = false;
              break;
            }
          }
          if (father_has_no_son_with_gene) { //if there is no gene in sons of the father, prune.
            if (not father->hasName())
              father->setName("Noname" + std::to_string(
                  noname_number++)); //if the father has no name we are going to have a leaf without name...

            if (verbose) {
              std::cout << "\t\tFather " << father->getName() << "(" << father << ")" << " has no descendance in genes"
                        << std::endl;
            }
            for (auto& son: sons) {
              // Then, delete nodes in the species_leafs list and in the tree
              for (auto it = species_leafs.begin(); it != species_leafs.end(); it++) {
                if ((*it) == son and current_species != (*it) and (*it) != father) {
                  species_leafs.erase(it);
                  break;
                }
              }
              speciestree.deleteNode(son);
              std::string name;
              if (verbose) {
                if (son->hasName())
                  name = son->getName();
                else if (son == speciestree.getRoot())
                  name = "the root";
                else
                  name = "";
              }
              if (verbose) std::cout << "\t\t\tSon: " << name << " (" << son << ")" << " erased." << std::endl;
            }
            // Then add the father to the leaf list because it is now a leaf !
            species_leafs.push_back(father);
          }
        } else { //if current_node is a leaf of the root
          //Reroot the tree.
          if(not genemap.hasGenes(speciestree.getRoot())) {
            speciestree.deleteNode(current_species);
            if(verbose) std::cout << "\t\t\tNode " << current_species << " erased." << std::endl;
            if (speciestree.getSons(father).size() == 1) {
              //if the root has no more than one son, reroot the tree to the son
              auto newroot = speciestree.getSons(father)[0];
              speciestree.rootAt(newroot);
              speciestree.deleteNode(father);
              if (verbose) std::cout << "\tTree is rerooted." << std::endl;
            }
          }
        }
      }
    }
    if(verbose) std::cout << "\t...leaf " << current_species << " done." << std::endl;
    species_leafs.erase(species_leafs.begin()); //remove the current_species of the vector.
  }

  // SECOND STEP pick leaves off the root
  bool end_of_root_pruning = false;
  while(not end_of_root_pruning){
    std::shared_ptr<bpp::PhyloNode> root = speciestree.getRoot();
    if(not genemap.hasGenes(root)) {
      auto sons = speciestree.getSons(root);
      end_of_root_pruning = true;
      for (auto son: sons) { //delete of leaves without gene.
        if (speciestree.isLeaf(son) and not genemap.hasGenes(son)) {
          end_of_root_pruning = false;
          speciestree.deleteNode(son);
        }
      }
      sons = speciestree.getSons(root);
      if (sons.size() == 1) {
        decltype(root) newroot = sons[0];
        speciestree.rootAt(newroot);
        speciestree.deleteNode(root);
        if (verbose) std::cout << "\tTree is rerooted." << std::endl;
      } else {
        end_of_root_pruning = true;
      }
    } else {
      end_of_root_pruning = true;
    }
  }
  if(verbose) std::cout << "Resulting tree: " << speciestree << std::endl;

  //resetNodeIdInPostOrder(speciestree);
}


std::shared_ptr<bpp::PhyloTree>
PhyloTreeToolBox::cloneTree(const bpp::PhyloTree &tree, const Node &root, const bool keep_branch_lengths, const bool verbose) {
  std::shared_ptr<bpp::PhyloTree> newTree(new bpp::PhyloTree());

  auto newTree_root = root != nullptr? root : tree.getRoot();

  std::unordered_set<Node> created_nodes;
  _recursiveTreeCloning_(tree, newTree_root, nullptr, *newTree, keep_branch_lengths, created_nodes, verbose);
  if(verbose) std::cout << "...done: " << *newTree << std::endl;

  //resetNodeIdInPostOrder(*newTree);

  newTree->rootAt(newTree_root);

  if(root == tree.getRoot() or root == nullptr) {
    assert(newTree->getAllEdges().size() == tree.getAllEdges().size() and
           newTree->getAllNodes().size() == tree.getAllNodes().size());
  }

  newTree->setName(tree.getName() + "_copy");

  assert(PhyloTreeToolBox::allEdgesAreCorrect(*newTree));
  return newTree;

}


std::vector<std::string> PhyloTreeToolBox::getLeavesNames(const bpp::PhyloTree &tree, const Node &root) {
  auto leaves = getLeaves(tree, root);
  return getNodeNames(leaves.begin(), leaves.end());
}

void PhyloTreeToolBox::_recursiveTreeCloning_(const bpp::PhyloTree &tree, const std::shared_ptr<bpp::PhyloNode> &current_node,
                                              const std::shared_ptr<bpp::PhyloNode> &father_node, bpp::PhyloTree &newtree,
                                              const bool keep_branch_lengths,
                                              std::unordered_set<Node>& created,
                                              const bool verbose) {
  if(verbose) std::cout << "\tVisiting " << current_node << "." << std::endl;
  if(not father_node) {
    if(verbose) std::cout << "\t\tCreating the root of the tree: " << current_node << "." << std::endl;
    /*
    if(newtree.isRooted()) {
      std::cerr << "Error " << __FUNCTION__ << ": tree is already rooted." << std::endl;
      assert(father_node != nullptr and not newtree.isRooted());
    }
     */
    assert(current_node);
    newtree.createNode(current_node);
    created.emplace(current_node);

    if(tree.hasIndex(current_node))
      newtree.setNodeIndex(current_node, tree.getNodeIndex(current_node));
  } else {
    assert(current_node);
    assert(father_node);
    assert(father_node);

    if(verbose) std::cout << "\t\tAdd node " << current_node << " son of " << father_node << "." << std::endl;

    // Create the new branch between father and son.
    std::shared_ptr<bpp::PhyloBranch> branch;

    if(not keep_branch_lengths) {
      branch = std::shared_ptr<bpp::PhyloBranch>(new bpp::PhyloBranch());
    } else {
      // Copy all properties of the original branch.
      if(not tree.getEdgeLinking(father_node, current_node)){
        std::cerr << "Error in tree named: " << tree.getName() << " rooted at " << tree.getRoot() << std::endl;
        std::cerr << tree << std::endl;
        std::cerr << "> Impossible to find edge between father " << father_node << " and son " << current_node << std::endl;
        std::cerr << "> With ids : " << tree.getNodeIndex(father_node) << " and " << tree.getNodeIndex(current_node) << std::endl;

        std::cerr << "=Summary of the tree===============================================================" << std::endl;

        std::cerr << "Nodes(" << tree.getAllNodes().size() << "): " << std::endl;
        auto tnodes = tree.getAllNodes();
        for(auto tnode: tnodes) std::cerr << tree.getNodeIndex(tnode) << ": " << tnode << std::endl;

        std::cerr << "Edges(" << tree.getAllEdges().size() << "): " << std::endl;
        auto edges = tree.getAllEdges();
        for(auto edge: edges) std::cerr << tree.getEdgeIndex(edge) << ": " << edge << std::endl;
        std::cerr << "=End summary=======================================================================" << std::endl;
        exit(EXIT_FAILURE);
      }
      branch = std::shared_ptr<bpp::PhyloBranch>(tree.getEdgeLinking(father_node, current_node)->clone());
    }

      // Then create node linking father and son.
    newtree.createNode(father_node, current_node, branch);
    assert(newtree.getEdgeToFather(current_node));
    created.emplace(current_node);

    if(tree.hasIndex(current_node)) {
      newtree.setNodeIndex(current_node, tree.getNodeIndex(current_node));
      newtree.setEdgeIndex(branch, tree.getNodeIndex(current_node));
    }
  }

  std::vector<Node> sons;
  for(auto son: tree.getOutgoingNeighbors(current_node)){
    if(created.find(son) == created.end())
      sons.push_back(son);
  }

  if(sons.size() > 0){
    for(auto son: sons) {
       _recursiveTreeCloning_(tree, son, current_node, newtree, keep_branch_lengths, created, verbose);
    }
  }
}

void PhyloTreeToolBox::sortNodesByIndex(std::vector<std::shared_ptr<bpp::PhyloNode>> &nodes, const bpp::PhyloTree &tree,
                                        const bool ascending) {
  std::sort(nodes.begin(), nodes.end(),
            [&tree, &ascending](const std::shared_ptr<bpp::PhyloNode>& node1, const std::shared_ptr<bpp::PhyloNode>& node2){
              return (ascending) == (tree.getNodeIndex(node1) < tree.getNodeIndex(node2)); //xnor operation
            }
  );
}

void PhyloTreeToolBox::sortNodesByNameLength(std::vector<std::shared_ptr<bpp::PhyloNode>> &nodes, const bool ascending) {
  std::sort(nodes.begin(), nodes.end(),
            [&ascending](const std::shared_ptr<bpp::PhyloNode>& node1, const std::shared_ptr<bpp::PhyloNode>& node2){
              std::string str1 = (node1->hasName()) ? node1->getName() : "";
              std::string str2 = (node2->hasName()) ? node2->getName() : "";
              return (ascending) == (str1.size() < str2.size()); //xnor operation
            }
  );
}

void PhyloTreeToolBox::addNodeWithSonsInPostOrder(bpp::PhyloTree &tree /// Tree to complete.
    , const std::shared_ptr<bpp::PhyloNode> &father /// Father node.
    , const std::shared_ptr<bpp::PhyloNode> &son_left /// Left son.
    , const std::shared_ptr<bpp::PhyloNode> &son_right /// Right son.
    , const std::shared_ptr<bpp::PhyloBranch> &son_left_branch /// Optional, PhyloBranch (in shared_ptr) between father and son left.
    , const std::shared_ptr<bpp::PhyloBranch> &son_right_branch /// Optional, PhyloBranch (in shared_ptr) between father and son right.
    , const bool verbose /// Print operations if true.
) {

  if(verbose) {
    std::cout << "\tAdd bifurcation: " << father << std::endl;
    if(son_left_branch)
      std::cout << "\t\t\t\t|_(" << (son_left_branch->hasLength() ? std::to_string(son_left_branch->getLength()) : "") << ")_ " << son_left << ((isArtificalGene(son_left)) ? "*" : "") << std::endl;
    else
      std::cout << "\t\t\t\t|__ " << son_left << ((isArtificalGene(son_left)) ? "*" : "") << std::endl;
    if(son_right_branch)
      std::cout << "\t\t\t\t|_(" << (son_right_branch->hasLength() ? std::to_string(son_right_branch->getLength()) : "") << ")_ " << son_right << ((isArtificalGene(son_right)) ? "*" : "") << std::endl;
    else
      std::cout << "\t\t\t\t|__ " << son_right << ((isArtificalGene(son_right)) ? "*" : "") << std::endl;
  }

  if(tree.hasFather(son_left)) {
    Node old_father = tree.getFather(son_left);
    tree.unlink(old_father, son_left);
  }

  if(tree.hasFather(son_right)) {
    Node old_father = tree.getFather(son_right);
    tree.unlink(old_father, son_right);
  }
  // get the current root
  Node root;
  if(not tree.isRooted())
    tree.rootAt(father);
  root = tree.getRoot();

  // link children to the father.
  if(son_left == root)
    tree.link(son_left, father, son_left_branch);
  else
    tree.link(father, son_left, son_left_branch);

  if(son_right == root)
    tree.link(son_right, father, son_right_branch);
  else
    tree.link(father, son_right, son_right_branch);

  // Change root for the most ancestral element (father).
  if(son_left == root or son_right == root) {
    tree.rootAt(father);
  }

  // Set properties of the father node
  bool isArtificial = isArtificalGene(son_left) and isArtificalGene(son_right);
  if(father->hasProperty("isArtificialGene")){
    auto& property = *dynamic_cast<IsArtificialGene*>(father->getProperty("isArtificialGene"));
    if(property.get() != isArtificial)
      property.set(isArtificial);
  } else {
    father->setProperty("isArtificialGene", IsArtificialGene(isArtificial));
  }
}

double
PhyloTreeToolBox::getDistanceBetweenTwoNodes(const Node &nodeA, const Node &nodeB,
                                             const bpp::PhyloTree &tree) {
  double distance = 0.0;
  auto path = tree.getEdgePathBetweenTwoNodes(nodeA, nodeB);
  for(auto& edge: path){
    if(edge->hasLength())
      distance += fabs(edge->getLength());
    else {
      std::cerr << __FUNCTION__ << " error (at " << __FILE__ << "): edge has no length." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  return distance;
}

bool
PhyloTreeToolBox::hasDistanceBetweenTwoNodes(const Node &nodeA, const Node &nodeB, const bpp::PhyloTree& tree) {
  auto path = tree.getEdgePathBetweenTwoNodes(nodeA, nodeB);
  for(auto& edge: path){
    if(not edge->hasLength()) {
      return false;
    }
  }
  return true;
}

std::shared_ptr<bpp::PhyloNode>
PhyloTreeToolBox::getNodeFromName(const bpp::PhyloTree &tree, const std::string& name, const bool case_sensitive, const bool leaves_only){
  std::vector<Node> nodes;
  if(leaves_only)
    nodes = tree.getAllLeaves();
  else
    nodes = tree.getAllNodes();
  return getNodeFromName(nodes, name, case_sensitive);
}

std::shared_ptr<bpp::PhyloNode>
PhyloTreeToolBox::getNodeFromName(const std::vector<Node> nodes, const std::string& name, const bool case_sensitive){
  for(auto it = nodes.begin(); it != nodes.end(); it++){
    // If the name is the same to an other in the tree, return its index.
    auto& node = *it;
    if(node->hasName()) {
      if (Utils::string_comp(node->getName(), name, case_sensitive))
        return node; // Node found
    } else {
      std::cerr << __FUNCTION__ << " error (at " << __FILE__ << " ): there is no name associated with the node." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  return nullptr; //FAIL
}

bool PhyloTreeToolBox::node_name_match_with_pattern(const std::string& node_name, const std::string& pattern, const bool case_sensitive) {
  bool exactName = true;
  bool name_extension_after = true;

  std::string PATTERN = pattern;
  // If there is a star somewhere (at the beginning or the end), its a pattern.
  if(pattern.front() == '*'){
    name_extension_after = false;
    PATTERN.erase(PATTERN.begin(), PATTERN.begin()+1);
    exactName = false;
  } else if(pattern.back() == '*') {
    PATTERN.erase(PATTERN.end() - 1, PATTERN.end());
    exactName = false;
  }

  if(node_name.size() >= (PATTERN.size())) {
    if(exactName){
      if(Utils::string_comp(PATTERN, node_name, case_sensitive)){
        return true;
      }
    } else {
      if (not name_extension_after) {
         return Utils::string_comp(PATTERN, node_name.substr(node_name.size() - PATTERN.size(), node_name.size()), case_sensitive);
      } else if (name_extension_after) {
        return Utils::string_comp(PATTERN, node_name.substr(0, PATTERN.size()), case_sensitive);
      }
    }
  }

  return false;
}

std::list<std::size_t>
PhyloTreeToolBox::getNodeIndexesFromNamePattern(const std::vector<Node>& nodes, const std::string& pattern, const bool case_sensitive){
  // Res contains gene candidates which matches with the std::string pattern
  std::list<std::size_t> res;

  for(auto it = nodes.begin(); it != nodes.end(); it++){
    const Node& node = *it;
    std::size_t i = (std::size_t) std::distance(nodes.begin(), it);
    /* If the name is the same to an other in the tree, return its index. */
    if(node->hasName()) {
      node_name_match_with_pattern(node->getName(), pattern, case_sensitive);
      res.push_back(i); // Add i.
    } else {
      std::cerr << __FUNCTION__ << " error: there is no name associated with the leaf." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  return res;
}

std::vector<std::shared_ptr<bpp::PhyloNode>>
PhyloTreeToolBox::getNodesFromNamePattern(const std::vector<Node>& nodes, const std::string& pattern, const bool case_sensitive){
  auto indexes_list = getNodeIndexesFromNamePattern(nodes, pattern, case_sensitive);
  std::vector<std::shared_ptr<bpp::PhyloNode>> res(indexes_list.size());
  auto indexes_list_iterator = indexes_list.begin();
  std::generate(res.begin(), res.end(), [&indexes_list_iterator, &nodes]{ return nodes.at(*(indexes_list_iterator++)); });

  assert(res.size() == indexes_list.size());

  return res;
}

Node PhyloTreeToolBox::getCommonAncestor(const bpp::PhyloTree &tree, const Node &nodeA, const Node &nodeB){
  std::list<Node> visited; //each node visited is added to the list.
  Node node = nodeA; //first node visited;
  assert(tree.isRooted());

  // First, visit each ancestor of nodeA
  while(node != tree.getRoot()) {
    visited.push_back(node);
    node = tree.getFather(node);
    if (node == nodeB) {
      return node;
    }
  }

  // Then, visit each ancestor of nodeB and check if the ancestor has been visited.
  node = nodeB;
  while(node != tree.getRoot()) {
    if(std::find(visited.begin(), visited.end(), node) != visited.end())
      return node;
    node = tree.getFather(node);
  }

  return node;
}

Node PhyloTreeToolBox::getLastCommonAncestor(std::vector<Node> nodes, const bpp::PhyloTree &tree) {
  Node& node = nodes.front();
  for(std::size_t i = 1; i < nodes.size(); i++){
    if(node == tree.getRoot()) {
      return node;
    }
    node = getCommonAncestor(tree, node, nodes.at(i));
  }
  return node;
}

std::list<Node>
PhyloTreeToolBox::getNodesInOrderTraversalRecursive_list(const bpp::PhyloTree &tree, const Node &root) {
  return InOrderTraversalRecursive_list(tree, [](const Node& node){ return true; }, root);
}

std::list<Node>
PhyloTreeToolBox::getNodesInPreOrderTraversalRecursive_list(const bpp::PhyloTree &tree, const Node &root) {
  return PreOrderTraversalRecursive_list(tree, [](const Node& node){ return true; }, root);
}

std::list<Node>
PhyloTreeToolBox::getNodesInPostOrderTraversalRecursive_list(const bpp::PhyloTree &tree, const Node &root) {
  return PostOrderTraversalRecursive_list(tree, [](const Node& node){ return true; }, root);
}

std::list<Node>
PhyloTreeToolBox::getNodesInPostOrderTraversalIterative_list(const bpp::PhyloTree &tree, const Node &root) {
  return getNodesInPostOrderTraversalIterative_list(tree, [](const Node &node) { return true; }, root);
}

std::vector<Node> PhyloTreeToolBox::getNodeBrothers(const bpp::PhyloTree &tree, const Node &node) {
  /// Node brothers are sharing the same father.
  std::vector<Node>  siblings = tree.getSons(tree.getFather(node));
  std::vector<Node> res;
  res.reserve(siblings.size() - 1);
  std::copy_if(siblings.begin(), siblings.end(), res.begin(),
               [&node](const Node& e) {
                 return e != node;
               });
  return res;
}

void PhyloTreeToolBox::resetNodeIdInPostOrder(bpp::PhyloTree &tree) {
  auto nodes = getNodesInPostOrderTraversalIterative_list(tree, tree.getRoot());
  std::size_t index = 0;
  for(auto node: nodes){
    tree.setNodeIndex(node, index);
    if(tree.hasFather(node)){
      tree.setEdgeIndex(tree.getEdgeToFather(node), index);
    }
    index++;
  }
}

std::list<Node> PhyloTreeToolBox::mergeWeakBranches(bpp::PhyloTree &tree, const double threshold, const bool inferior_only) {
  // Define a functor which indicates when a branch (and the node) is weak.
  std::list<Node> deletedNodes;
  auto isweak = [&threshold, &deletedNodes, &inferior_only](const Node& node, const bpp::PhyloTree& tree) {
    if(tree.getRoot() != node and (not tree.isLeaf(node))) {
      auto edge = tree.getEdgeToFather(node);
      if(edge->hasBootstrapValue()){
        auto bootstrap = edge->getBootstrapValue();
        if(Utils::double_equal_or_inferior(bootstrap, threshold)) {
          if(inferior_only
             and Utils::double_equivalence(bootstrap, threshold)) {
            return false;
          } else {
            deletedNodes.push_back(node);
            return true;
          }
        }
        else {
          //edge->removeProperty("bootstrap");
          return false;
        }
      } else {
        // Case where the is no support associated with the branch
        if(CONTRACT_BRANCHES_WITHOUT_SUPPORT) {
          if(Utils::double_equal_or_inferior(DEFAULT_BRANCH_SUPPORT_VALUE, threshold)) {
            if(inferior_only
               and Utils::double_equivalence(DEFAULT_BRANCH_SUPPORT_VALUE, threshold)) {
              return false;
            } else {
              deletedNodes.push_back(node);
              return true;
            }
          }
        }
      }
    }
    return false;
  };

  removeNodesInPostOrderTraversal(tree, tree.getRoot(), isweak, false);

  PhyloTreeToolBox::resetNodeIdInPostOrder(tree);

  return deletedNodes;
}

void PhyloTreeToolBox::removeNode(bpp::PhyloTree &tree, const Node& node, const bool verbose) {
  Node current = node;
  if(not tree.isLeaf(current) and tree.getRoot() != node){
    // If the current node is internal, link all sons to their grandpa.
    auto grandpa = tree.getFather(current);
    auto sons = tree.getSons(current);
    bool set_branch_length = false;
    for (auto &son: sons) {
      auto sonEdgeToFather = tree.getEdgeToFather(son);
      auto grandpa_to_current_son_path = tree.getEdgePathBetweenTwoNodes(grandpa, son);
      for(auto& edge_path: grandpa_to_current_son_path)
        if(edge_path->hasLength()){
          set_branch_length = true;
          break;
        }

      std::shared_ptr<bpp::PhyloBranch> branch = nullptr;
      if(set_branch_length)
        branch = std::shared_ptr<bpp::PhyloBranch>(
          new bpp::PhyloBranch(getDistanceBetweenTwoNodes(grandpa, son, tree)));

      //std::cout << "Change for: " << grandpa << " --(" << branch->getLength() << ")-- " << son << std::endl;
      tree.link(grandpa, son, branch);
      if(tree.hasIndex(sonEdgeToFather))
        tree.setEdgeIndex(branch, tree.getEdgeIndex(sonEdgeToFather));
    }
  }
  tree.deleteNode(current);
  if(verbose) std::cout << __FUNCTION__ << ": node " << current << " deleted." << std::endl;
}

/*
void PhyloTreeToolBox::pruneTree_v1(bpp::PhyloTree &speciestree,
                                    const SpeciesGeneMap &genemap,
                                    const bool verbose) {
  struct isNotAssociated {
    SpeciesGeneMap genemap_;
    PhyloTreeToolBox toolbox_;
    isNotAssociated(const SpeciesGeneMap& genemap): genemap_(genemap) {};
    bool operator()(const Node& node, const bpp::PhyloTree& tree) {
      bool result = true;
      if(tree.isLeaf(node)){
        if(not genemap_.hasGenes(node)){
          if(tree.getRoot() != tree.getFather(node)) {
            auto brothers = toolbox_.getNodeBrothers(tree, node);
            for (auto &brother: brothers) {
              if (genemap_.hasGenes(brother)) {
                result = false;
                break;
              }
            }
          } else {
            result = false; // All leaves connected with the root and not associated with a gene is deleted.
          }
        } else {
          result = false;
        }
      } else {
        result = false;
      }
      return result;
    }
  };
  removeNodesInPostOrderTraversal(speciestree, speciestree.getRoot(), isNotAssociated(genemap), false);
}
*/

void PhyloTreeToolBox::removeArtificialGenes(bpp::PhyloTree &tree) {
  applyInPOT(tree,
             [&tree](const Node &node) {
               if (not tree.isLeaf(node)) {
                 if (tree.getSons(node).size() == 1) {
                   if(tree.getRoot() == node){
                     tree.rootAt(tree.getSons(node).at(0));
                   }
                   PhyloTreeToolBox::removeNode(tree, node);
                 }
               } else {
                 if(isArtificalGene(node)){
                   PhyloTreeToolBox::removeNode(tree, node);
                 }
               }
             });
}


bool PhyloTreeToolBox::isBinary(const bpp::PhyloTree &tree, const Node &root) {

  return not hasNode(tree, [&tree](const Node& node){
    if(not tree.isLeaf(node)){
      std::vector<Node> sons = tree.getSons(node);
      return sons.size() != 2;
    }
    return false;
  });
}

bool PhyloTreeToolBox::compare_trees_with_days_algorithm(const bpp::PhyloTree &treeA, const bpp::PhyloTree &treeB,
                                                         const bool verbose) {
  /// @note See Thomas G. Kristensen, Study of a Simple Pruning Strategy with Days Algorithm.
  auto nodesA = getNodesInPostOrderTraversalIterative_list(treeA);
  auto nodesB = getNodesInPostOrderTraversalIterative_list(treeB);
  auto nodesB_leaves = getNodesInPostOrderTraversalIterative_list(treeB, [&treeB](const Node &node) {
    return not treeB.isLeaf(node);
  });

  if(nodesA.size() != nodesB.size()) {
    if(verbose) std::cout << __FUNCTION__ << ": the two trees has different number of nodes (" << nodesA.size() << " and " << nodesB.size() << ")." << std::endl;
    return false;
  }

  std::map<std::string, std::size_t> labels; // the key is the node name, and the value an index given by the algorithm.

  // each internal node is associated with a pair of indexes which is an interval describing leaves under.
  // these intervals are compared to the other in the second tree.
  std::unordered_map<Node, std::pair<std::size_t, std::size_t>> intervalsA;
  std::unordered_map<Node, std::pair<std::size_t, std::size_t>> intervalsB;

  std::size_t i = 0; // will give the index given by the algorithm.
  for(auto& node: nodesA) {
    if (treeA.isLeaf(node))
      labels[node->getName()] = i++;
    else {
      //Because of Post-Order traversal, we can compute intervals in the same loop.
      auto sons = treeA.getSons(node);
      if (sons.size() != 2) std::cout << "First tree is not valid." << std::endl;
      assert(sons.size() == 2);
      auto son_left = sons.at(0);
      auto son_right = sons.at(1);

      std::size_t min_son_left;
      std::size_t max_son_left;
      std::size_t min_son_right;
      std::size_t max_son_right;

      if (treeA.isLeaf(son_left)) {
        min_son_left = labels.at(son_left->getName());
        max_son_left = labels.at(son_left->getName());
      } else {
        auto son_interval = intervalsA.at(son_left);
        min_son_left = son_interval.first;
        max_son_left = son_interval.second;
        if (min_son_left > max_son_left)
          std::swap(min_son_left, max_son_left);
      }
      if (treeA.isLeaf(son_right)) {
        min_son_right = labels.at(son_right->getName());
        max_son_right = labels.at(son_right->getName());
      } else {
        auto son_interval = intervalsA.at(son_right);
        min_son_right = son_interval.first;
        max_son_right = son_interval.second;
        if (min_son_right > max_son_right)
          std::swap(min_son_right, max_son_right);
      }

      std::size_t left = min_son_left;
      std::size_t right = max_son_right;
      if(max_son_left < min_son_right){
        left = min_son_left;
        right = max_son_right;
      }
      else if(max_son_right < min_son_left){
        left = min_son_right;
        right = max_son_left;
      }
      else if(min_son_left < min_son_right and max_son_left < max_son_right and max_son_left > min_son_right){
        left = min_son_left;
        right = max_son_right;
      }
      else if(min_son_right < min_son_left and max_son_right < max_son_left and max_son_right > min_son_left){
        left = min_son_right;
        right = max_son_left;
      }

      if (left < right)
        intervalsA[node] = std::pair<std::size_t, std::size_t>(left, right);
      else
        intervalsA[node] = std::pair<std::size_t, std::size_t>(right, left);
    }
  }

  for(auto& node: nodesB_leaves) {
    //Because of Post-Order traversal, we can compute intervals in the same loop.
      auto sons = treeB.getSons(node);
      if(sons.size() != 2) std::cout << "Second tree is not valid." << std::endl;
      assert(sons.size() == 2);
      auto son_left = sons.at(0);
      auto son_right = sons.at(1);
      std::size_t min_son_left;
      std::size_t max_son_left;
      std::size_t min_son_right;
      std::size_t max_son_right;
      if(treeB.isLeaf(son_left)){
        min_son_left = labels.at(son_left->getName());
        max_son_left = labels.at(son_left->getName());
      }else{
        auto son_interval = intervalsB.at(son_left);
        min_son_left = son_interval.first;
        max_son_left = son_interval.second;
        if(min_son_left  > max_son_left)
          std::swap(min_son_left, max_son_left);
      }
      if(treeB.isLeaf(son_right)){
        min_son_right = labels.at(son_right->getName());
        max_son_right = labels.at(son_right->getName());
      }else{
        auto son_interval = intervalsB.at(son_right);
        min_son_right = son_interval.first;
        max_son_right = son_interval.second;
        if(min_son_right > max_son_right)
          std::swap(min_son_right, max_son_right);
      }

      std::size_t left = min_son_left;
      std::size_t right = max_son_right;
      if(max_son_left < min_son_right){
        left = min_son_left;
        right = max_son_right;
      }
      else if(max_son_right < min_son_left){
        left = min_son_right;
        right = max_son_left;
      }
      else if(min_son_left < min_son_right and max_son_left < max_son_right and max_son_left > min_son_right){
        left = min_son_left;
        right = max_son_right;
      }
      else if(min_son_right < min_son_left and max_son_right < max_son_left and max_son_right > min_son_left){
        left = min_son_right;
        right = max_son_left;
      }

      if(left < right)
        intervalsB[node] = std::pair<std::size_t, std::size_t>(left, right);
      else
        intervalsB[node] = std::pair<std::size_t, std::size_t>(right, left);

      auto interB = intervalsB.at(node);
      auto it_interA = std::find_if(std::begin(intervalsA), std::end(intervalsA),
                                    [&interB](const std::pair<Node, std::pair<std::size_t, std::size_t>>& pairA){
                                      auto interA = pairA.second;
                                      if(interA.first == interB.first and interA.second == interB.second){
                                        return true;
                                      }
                                      return (interA.first == interB.second and interA.second == interB.first);
                                    });
      if(it_interA == std::end(intervalsA)){
        std::cout << "Intervals A: " << std::endl << intervalsA << std::endl;
        std::cout << "Intervals B: " << std::endl << intervalsB << std::endl;
        std::cout << "Labels: " << std::endl << labels << std::endl;
        std::cout << "Difference in ancestral node " << interB << "." << std::endl;
        std::cout << "Tree A: " << treeA << std::endl;
        std::cout << "Tree B: " << treeB << std::endl;
        return false;
      }
      intervalsA.erase(it_interA);
  }
  return intervalsA.empty();
}

bool PhyloTreeToolBox::hasArtificialGenes(const bpp::PhyloTree &tree, const Node& root) {
  return hasNode(tree, [](const Node &node) { return PhyloTreeToolBox::isArtificalGene(node); }, root);
}

bool PhyloTreeToolBox::isArtificalGene(const Node &node) {
  if(node->hasProperty("isArtificialGene")){
    return dynamic_cast<IsArtificialGene*>(node->getProperty("isArtificialGene"))->get();
  }
  return false;
}

std::size_t PhyloTreeToolBox::getNodeOccurence(const Node &node) {
  if(node->hasProperty("occurence")){
    return dynamic_cast<NodeOccurence*>(node->getProperty("occurence"))->get();
  }
  std::cerr << "Error, there is no occurence of this gene." << std::endl;
  exit(EXIT_FAILURE);
}

std::list<Node> PhyloTreeToolBox::find_nodes_according_to_its_outdegree(const bpp::PhyloTree &tree,
                                                                        const std::size_t &outdegree, const Node &root
) {
  /// The third argument (const Node root) is
  /// optional, it is the root of the sub-tree in the tree, defaults is the root of the entire tree.
  return getNodesInPostOrderTraversalIterative_list(tree, [&tree, &outdegree](const Node &node) {
    if (not tree.isLeaf(node)) {
      if (tree.getSons(node).size() == outdegree) return true;
    }
    return false;
  }, root);
}

std::list<Node> PhyloTreeToolBox::findPolytomies(const bpp::PhyloTree &tree, const Node &root) {
  return getNodesInPostOrderTraversalIterative_list(tree, [&tree](const Node &node) {
    if (not tree.isLeaf(node)) {
      if (tree.getSons(node).size() > 2) return true;
    }
    return false;
  }, root);
}

std::list<Node> PhyloTreeToolBox::findBifurcations(const bpp::PhyloTree &tree, const Node& root) {
  return find_nodes_according_to_its_outdegree(tree, 2, root);
}

std::list<Node> PhyloTreeToolBox::findMonofurcations(const bpp::PhyloTree &tree, const Node& root) {
  return find_nodes_according_to_its_outdegree(tree, 1, root);
}

bool PhyloTreeToolBox::hasMonofurcations(const bpp::PhyloTree &tree, const Node& root) {
  return hasNode(tree,
                 [&tree](const Node &node) {
                   if (not tree.isLeaf(node))
                     return tree.getSons(node).size() == 1;
                   else return false;
                 },
                 root
  );
}

std::pair<std::list<Node>, std::list<Node>> PhyloTreeToolBox::getFurcations(const bpp::PhyloTree &tree, const Node& root) {
  std::list<Node> bifurcationRoots;
  std::list<Node> multifurcationRoots;
  //auto nodes = tree.getAllInnerNodes();
  //sortNodesByIndex(nodes, tree);
  auto nodes = getNodesInPostOrderTraversalIterative_list(tree, root);
  for (auto &node: nodes) {
    // If the node has more than 2 sons, it is a polytomy.
    if(not tree.isLeaf(node)) {
      if (tree.getSons(node).size() == 2) {
        bifurcationRoots.push_back(node);
      } else if (tree.getSons(node).size() > 2) {
        multifurcationRoots.push_back(node);
      }
    }
  }
  return std::pair<decltype(bifurcationRoots), decltype(multifurcationRoots)>(bifurcationRoots, multifurcationRoots);
}

std::map<Event, std::size_t>
PhyloTreeToolBox::getGeneEventsInBifurcation(const Node &node, const bpp::PhyloTree &genetree,
                                             const bpp::PhyloTree &speciestree, const SpeciesGeneMap &map) {
  assert(not genetree.isLeaf(node));

  std::map<Event, std::size_t> events;
  events[duplication] = 0;
  events[loss] = 0;

  // Get children of the current_gene
  auto node_children = genetree.getSons(node);
  assert(node_children.size() == 2);

  // Get species associated with the current gene (node)
  auto node_species = map.getAssociatedSpecies(node);

  // Get children species
  auto& left_child = node_children.at(0);
  auto& right_child = node_children.at(1);
  auto left_child_species = map.getAssociatedSpecies(left_child);
  auto right_child_species = map.getAssociatedSpecies(right_child);

  // Get species sons of the node species.
  auto species_children = speciestree.getSons(node_species);

  if(not genetree.isLeaf(node) and (node_species == left_child_species or node_species == right_child_species))
    events[duplication] += 1;

  for(auto& node_child: node_children) {
    // for leaf left and leaf right

    // get species associate with the leaf
    auto node_child_species = map.getAssociatedSpecies(node_child);

    // if there is no duplication recorded
    if(events.at(duplication) == 0){
      // while the species leaf do not correspond with a species son of the bifurcation root
      // add a loss and do it until we are in a species son of the bifurcation root.
      while(std::find(species_children.begin(), species_children.end(), node_child_species) == species_children.end()) {
        events[loss] += 1;
        if(node_child_species == speciestree.getRoot()) break;
        node_child_species = speciestree.getFather(node_child_species);
      }
    }
    else if(events.at(duplication) > 0){
      while(node_child_species != node_species) {
        events[loss] += 1;
        //assert(node_child_species != speciestree.getRoot());
        if(node_child_species == speciestree.getRoot()) break;
        node_child_species = speciestree.getFather(node_child_species);
      }
    }
  }

  return events;
}

std::list<Node> PhyloTreeToolBox::getLeaves(const bpp::PhyloTree &tree, const Node &root) {
  return getNodesInPostOrderTraversalIterative_list(tree,
                                                    [&tree](const Node& node) {
                                                      return tree.isLeaf(node);
                                                    },
                                                    root);
}

std::vector<std::string> PhyloTreeToolBox::getInternalNodeNames(const bpp::PhyloTree &tree, const Node &root) {
  auto i_nodes = getInternalNodes(tree, root);
  return getNodeNames(i_nodes.begin(), i_nodes.end());
}

std::list<Node> PhyloTreeToolBox::getInternalNodes(const bpp::PhyloTree &tree, const Node &root) {
  return getNodesInPostOrderTraversalIterative_list(tree,
                                                    [&tree](const Node& node) {
                                                      return not tree.isLeaf(node);
                                                    },
                                                    root);
}

Event PhyloTreeToolBox::getEvent(const bpp::PhyloTree &genetree, const SpeciesGeneMap &map, const Node &node) {
  if(genetree.isLeaf(node)){
    if(PhyloTreeToolBox::isArtificalGene(node)) {
      return Event::loss;
    } else {
      return Event::extant;
    }
  } else {
    auto sons = genetree.getSons(node);
    auto& son_left = sons.front();
    auto& son_right = sons.back();

    if(son_left == son_right){
      if(sons.size() == 1)
        std::cerr << "Error (PhyloTreeToolBox::getEvent): an internal node " << node << " has only one son." << std::endl;
      else
        std::cerr << "Error (PhyloTreeToolBox::getEvent): an internal node " << node << "  has repeated children." << std::endl;
    }

    auto son_left_species = map.getAssociatedSpecies(son_left);
    auto son_right_species = map.getAssociatedSpecies(son_right);

    if(son_left_species == son_right_species){
      return Event::duplication;
    } else {
      /*
      if(getEvent(genetree, map, son_left) == Event::loss or getEvent(genetree, map, son_right) == Event::loss){
        //return Event::speciationLoss;
        //} else
        */
      return Event::speciation;
    }
  }

  return Event::none;
}

bool PhyloTreeToolBox::allEdgesAreCorrect(const bpp::PhyloTree &tree) {
  return not PhyloTreeToolBox::hasNode(tree,
                                       [&tree](const Node& node){
                                         if(tree.getRoot() != node) {
                                           return not tree.getEdgeToFather(node);
                                         }
                                         return false; //if node is the root.
                                       }
  );
}


