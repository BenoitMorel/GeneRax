// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 28/12/16.
//

#include "IO.h"

// Include Bpp
#include <Bpp/BppString.h>

// Include Treerecs-code
#include <ale/tools/PhyloTreeToolBox.h>
#include <ale/tools/Utils.h>

/// bpp::Newick module.
const Newick IO::newick_ = Newick();
/// bpp::NHX module.
const Nhx IO::nhx_ = Nhx(false);

std::shared_ptr<bpp::PhyloTree> IO::nhxToPhyloTree(const std::string &description){
  return std::shared_ptr<bpp::PhyloTree>(nhx_.parenthesisToPhyloTree(description));
}

std::shared_ptr<bpp::PhyloTree>
IO::newickToPhyloTree(const std::string &description, bool bootstrap, const std::string &propertyName, bool withId,
                      bool verbose){
  return std::shared_ptr<bpp::PhyloTree>(newick_.parenthesisToPhyloTree(description, bootstrap, propertyName, withId, verbose));
}

std::shared_ptr<bpp::PhyloTree>
IO::readTreeFile(const std::string &filename, const Text_format format, const bool support, const bool check_branch_length,
                 const bool verbose){
  return IO::readTreesFile(filename, format, 0, support, check_branch_length, 0, verbose).front();
}

std::vector<std::shared_ptr<bpp::PhyloTree>>
IO::readTreesFile(const std::string &filename, const Text_format format, const int index, const bool support,
                  const bool check_branch_length, const bool printProgression, const bool verbose) {

  if(verbose) std::cout << "Read " << filename << " with " << format << " format." << std::endl;
  std::vector<std::shared_ptr<bpp::PhyloTree>> trees;


  // Open the file
  std::ifstream genefile(filename.c_str(), std::ios::in);
  try {
    if (not genefile.is_open()) {
      throw std::invalid_argument(filename + std::string(" does not exist."));
    }
  } catch(std::exception const& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  // Get file content and strip line feeds.
  std::string file_content((std::istreambuf_iterator<char>(genefile)), std::istreambuf_iterator<char>() ); // content of the file

  /*
  unsigned int corresponding_nchar_with_linefeed = 0;
  unsigned int nchar_linefeed = strlen(LINE_FEED);
  for( auto str_it = file_content.begin() ; str_it != file_content.end() ; str_it++) {
      if(corresponding_nchar_with_linefeed == nchar_linefeed) {
          assert(*(str_it -corresponding_nchar_with_linefeed) == LINE_FEED[0]);
          file_content.erase(str_it - corresponding_nchar_with_linefeed, str_it);
          corresponding_nchar_with_linefeed = 0;
      }

      assert(corresponding_nchar_with_linefeed < nchar_linefeed);

      if((*str_it) == LINE_FEED[corresponding_nchar_with_linefeed]){
          corresponding_nchar_with_linefeed++;
      } else {
          if(corresponding_nchar_with_linefeed > 0) {
              corresponding_nchar_with_linefeed = 0;
          }
      }
  }

  if(corresponding_nchar_with_linefeed == nchar_linefeed) {
      assert(*(file_content.end() -corresponding_nchar_with_linefeed) == LINE_FEED[0]);
      file_content.erase(file_content.end() - corresponding_nchar_with_linefeed, file_content.end());
      corresponding_nchar_with_linefeed = 0;
  }
  */

  file_content.erase(
      std::remove(file_content.begin(), file_content.end(), '\n')
      , file_content.end());

  #if defined _WIN32 || defined __CYGWIN__
  file_content.erase(
      std::remove(file_content.begin(), file_content.end(), '\r')
      , file_content.end());
  #endif

  Utils::trim_str(file_content);

  genefile.close();

  std::size_t number_of_trees = 0;
  if(format == Text_format::newick or format == Text_format::nhx) {
    number_of_trees = Utils::count(file_content, ";");
  } else if(format == Text_format::phyloxml) {
    number_of_trees = Utils::count(file_content, "phylogeny")/2;
  } else {
    std::cerr << "Error: not supported tree format." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(verbose) std::cout << "Number of trees found in " << filename << ": " << number_of_trees << std::endl;

  if(number_of_trees == 0) {
    throw std::runtime_error("IO: no tree read in " + filename + ".");
    return {};
  };

  if(index == -1) {
    trees.reserve(number_of_trees);
  }
  else if (index >= 0) {
    if(static_cast<std::size_t>(index) >= number_of_trees) {
      std::cerr << "Error: index " << index + 1
                << " is greater than the number of trees in " << filename
                << " (" << number_of_trees << ")." << std::endl;
      exit(EXIT_FAILURE);
    }
    trees.reserve(1);
  }

  if(format == Text_format::newick) {
    auto contents = Utils::splitString(file_content, ";");
    int i = 0;
    for(auto content: contents){
      if(index == -1 or (index >= 0 and i == index)) {
        auto tree = std::shared_ptr<bpp::PhyloTree>(IO::newickToPhyloTree(content + ";", support));
        if(tree) {
          trees.emplace_back(tree);
          if (verbose) std::cout << *trees.back() << std::endl;
        } else {
            std::cout << std::endl << "Error: bad format for \"" << content << "\"" << std::endl;
          throw std::runtime_error("IO: bad format in " + filename + ".");
        }
      }
      i++;
    }
  }
  else if(format == Text_format::nhx) {
    auto contents = Utils::splitString(file_content, ";");
    int i = 0;
    for(auto content: contents){
      if(index == -1 or (index >= 0 and i == index)) {
        auto tree = std::shared_ptr<bpp::PhyloTree>(IO::nhxToPhyloTree(content + ";"));
        if(tree) {
          trees.emplace_back(tree);
        }else{
          throw std::runtime_error("IO: bad format in " + filename + ".");
        }
      }
      i++;
    }
  }

  if(trees.size() == 0) {
    std::cerr << "No tree found in " << filename << " with " << format << " format." << std::endl;
    throw std::runtime_error("IO: no tree read in " + filename + ".");
    return {nullptr};
  };

  for(std::size_t i = 0; i < trees.size(); ++i) {
    // Read the content and control

    auto& tree = trees.at(i);

    assert(tree != nullptr);
    if(tree->getName().empty()){
      if(index >= 0)
        tree->setName("Tree" + std::to_string(index + 1));
      else
        tree->setName("Tree" + std::to_string(i + 1));
    }

    // Set a name to each internal node (only leaves has a name yet).
    preProcesses(*tree, check_branch_length, false);
  }

  return trees;
}


void IO::preProcesses(bpp::PhyloTree &tree
    , const bool check_branch_length /// Check if branch lengths are correct, if there is no branch length it will be replaced with a default value and warning will be printed.
    , const bool force_branch_support /// Add support when there is no support/bootstrap on a branch.
){
  // First of all: find if a leaf name is a digit or a number, because of each node will be named as the index defined
  // by Bio++...
  int increment = 0;
  for(auto node: tree.getAllLeaves()){
    node->setProperty("isArtificialGene", IsArtificialGene(false));
    if(Utils::stringIsNumber(node->getName())) {
      // std::cerr << "The tree can't have an elements named as a number: " << node->getName() << "." << std::endl;
      int num = std::stoi(node->getName());
      if(num > increment) {
        increment = (int) (num + 1);
      }
    }
  }

  // Then, set a name to the internal nodes, the name is going to be an integer which increments as the post order traversal.
  for(auto node: PhyloTreeToolBox::getInternalNodes(tree)){
    node->setProperty("isArtificialGene", IsArtificialGene(false));
    if (!node->hasName()) {
      //auto nodeIndex = tree.getNodeIndex(node);
      //node->setName(std::to_string(nodeIndex + increment));
      node->setName(std::to_string(increment));
      increment++;
    }

    if (force_branch_support and (node != tree.getRoot())) {
      auto edgeToFather = tree.getEdgeToFather(node);
      if (not edgeToFather->hasBootstrapValue()) {
        edgeToFather->setProperty(
            "bootstrap",
            bpp::Number<double>(DEFAULT_BRANCH_SUPPORT_VALUE));
      }
    }
  }

  if(check_branch_length) {
    bool tree_has_edge_without_length = false;
    std::size_t n_edges_without_length = 0;
    double init_edge_length_value = MINIMAL_BRANCH_LENGTH;
    for (auto edge: tree.getAllEdges()) {
      if (not edge->hasLength()) {
        if (not tree_has_edge_without_length) tree_has_edge_without_length = true;
        edge->setLength(init_edge_length_value);
        n_edges_without_length += 1;
      }
    }

    if (tree_has_edge_without_length) {
      std::cerr << "Warning: tree has " << n_edges_without_length << " branch"
                << (n_edges_without_length > 1 ? "s" : "")
                << " without length." << std::endl;
      std::cerr << "Branch length set to " << init_edge_length_value << std::endl;
    }
  }
}

std::string IO::PhyloTreeToNewick(const bpp::PhyloTree &tree){
  return newick_.treeToParenthesis(tree, false);
}

std::string IO::PhyloTreeToNewick(const bpp::PhyloTree &tree, bool bootstrap, const std::string &propertyName){
  return newick_.treeToParenthesis(tree, bootstrap, propertyName);
}

std::string IO::PhyloTreeToNhx(const bpp::PhyloTree &tree){
  return nhx_.treeToParenthesis(tree);
}


void IO::write(
    const bpp::PhyloTree &tree
    , std::ostream &os
    , const Text_format output
    , const std::string& description
    , const bool reconciliation_mode
){
  if(not description.empty() and (output != Text_format::phyloxml and output != Text_format::recphyloxml)) {
    os << "> " << description << std::endl;
  } else if(not description.empty()) {
    os << "<!-- " << description << " -->" << std::endl;
  }

  if(output == Text_format::nhx or output == Text_format::recnhx)
    nhx_.write(tree, os);
  else
    newick_.write(tree, os, false);
}

bool IO::exists(const std::string &filename) {
  return ( access( filename.c_str(), F_OK ) != -1 );
    /*!
  if (FILE *file = fopen(filename.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
     */
}

unsigned long IO::nlines(std::ifstream &is) {
  // Reset the file stream and start the mapping.
  is.clear();
  is.seekg(0, is.beg);

  // Compute the number of lines.
  unsigned long res;
  std::string line_content;
  for (res = 0; std::getline(is, line_content); ++res)
    ;

  // Reset the file stream and start the mapping.
  is.clear();
  is.seekg(0, is.beg);

  return res;
}
