// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 28/12/16.
//

#ifndef PHYLASOLVER_IO_H
#define PHYLASOLVER_IO_H

// Include libs
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <unordered_map>
#include <list>
#include <memory>
#include <stdexcept>

// Include Bpp
#include <Bpp/Phyl/Tree/PhyloTree.h>

#include <ale/Defines.h>

// Include containers
#include <ale/containers/NodeProperty.h>

//Include tools
#include "Nhx.h"
#include "Newick.h"

/// Format of tree's inputs/outputs.
enum class Text_format {
  newick,
  nhx,
  recnhx,
  phyloxml,
  recphyloxml,
  svg,
  unknown
};

inline Text_format operator++( Text_format& x ) {
  return x = static_cast<Text_format>((
      static_cast<int>(x) + 1));
}

/// Write an Text_format into an output stream.
inline std::ostream& operator<<(std::ostream& os, const Text_format& output_format) {
  switch (output_format) {
    case(Text_format::nhx): {
      os << "NHX";
      break;
    }
    case(Text_format::recnhx): {
      os << "NHX+";
      break;
    }
    case(Text_format::newick): {
      os << "Newick";
      break;
    }
    default:
      os << "Unknown";
  }
  return os;
}

/*!
 * @class IO
 * @brief IO provides functions to load a bpp::PhyloTree or write it into a Newick, Nhxfile.
 * @details Trees are returned as std::shared_ptr<bpp::PhyloTree> and each node is named.
 */
class IO {

private:

protected:
  /// Newick module to write/ read Newick parenthesis format.
  static const Newick newick_;
  /// NHX module to write/ read NHX parenthesis format.
  static const Nhx nhx_;

  /// \brief Set a name for each unnamed internal node, check if the tree has branch lengths in options.
  /// \return Nothing but edit the given tree.
  static void preProcesses(bpp::PhyloTree& tree /// Tree to modify.
      , const bool check_branch_length = true /// Check if tree has correct branch lengths (default = true).
      , const bool force_branch_support = false /// Set branch support when the branch has no one.
  );

  /// \brief Rename homonyms in tree by adding a number.
  /// \tparam NodesIterator Iterator on std::shared_ptr<bpp::PhyloNode>
  /// \param nodes_begin
  /// \param nodes_end
  template<typename NodesIterator>
  static void renameNodeHomonyms(
      const NodesIterator &nodes_begin
      , const NodesIterator &nodes_end
  ) {
    std::unordered_map<std::string, std::list<Node>> nodes_by_name;

    // First we classified nodes by their names
    NodesIterator nodes_it;
    for(nodes_it = nodes_begin ; nodes_it != nodes_end ; nodes_it++) {
      auto& node = *nodes_it;
      if(node->hasName()) {
        std::string gene_name = node->getName();
        nodes_by_name[gene_name].push_back(node);
      }
    }

    // Then, if one list has a length > 1, we rename all genes in this list.
    auto nbn_it = nodes_by_name.begin();
    for(; nbn_it != nodes_by_name.end(); nbn_it++) {
      if(nbn_it->second.size() > 1) {
        unsigned int occurrence = 0;
        for (auto node: nbn_it->second) {
          auto node_name = node->getName();
          node_name += ("_" + std::to_string(occurrence));
          node->setName(node_name);
          occurrence++;
        }
      }
    }
  }

public:
/****************
 * LOAD
 */
  /// Check file existence.
  static bool exists(const std::string& filename);

  /// Compute number of lines in a file. The file is going to be read and reinit.
  static unsigned long nlines(std::ifstream& is);

  /// Translates newick std::string to bpp::PhyloTree.
  static std::shared_ptr<bpp::PhyloTree> nhxToPhyloTree(const std::string& description);

  /// Translates newick std::string to bpp::PhyloTree.
  static std::shared_ptr<bpp::PhyloTree> newickToPhyloTree(const std::string& description, bool bootstrap = false, const std::string& propertyName = "", bool withId = false, bool verbose = false);

  /// Read a tree file and return a std::shared_ptr<bpp::PhyloTree>.
  static std::shared_ptr<bpp::PhyloTree> readTreeFile(const std::string &filename /// name of the newick file.
      , const Text_format format = Text_format::newick/// file in nhx format (default = false).
      , const bool support = false /// get support values.
      , const bool check_branch_length = true /// check branch lengths. Prints warning when there is no branch length.
      , const bool verbose = false /// verbose mode.
  );

  /// Read a Newick file and return a vector of std::shared_ptr<bpp::PhyloTree>.
  static std::vector<std::shared_ptr<bpp::PhyloTree>>
  readTreesFile(const std::string &filename
                , const Text_format format = Text_format::newick
                , const int index = -1
                , const bool support = true
                , const bool check_branch_length = true
                , const bool printProgression = true
                , const bool verbose = false
  );


/****************
 * SAVE
 */
  /// \brief Translate bpp::PhyloTree into newick's parenthesis format.
  /// \return Std::string which contains the tree in the newick parenthesis format.
  static std::string PhyloTreeToNewick(const bpp::PhyloTree &tree);

  /// \brief Translate bpp::PhyloTree into newick's parenthesis format.
  /// @param tree The tree to convert.
  /// @param bootstrap Tell is bootstrap values must be writen.
  ///   If so, the content of the property with name "bootstrap" will be written as bootstrap value.
  ///   The property should be a Number<double> object.
  ///   Otherwise, the content of the property with name 'propertyName' will be written.
  ///   In this later case, the property should be a String object.
  /// @param propertyName The name of the property to use. Only used if bootstrap = false.
  /// \return Std::string which contains the tree in the newick parenthesis format.
  static std::string PhyloTreeToNewick(const bpp::PhyloTree& tree, bool bootstrap, const std::string& propertyName);

  /// \brief Translate bpp::PhyloTree into nhx parenthesis format.
  /// \return Std::string which contains the tree in the NHX parenthesis format.
  static std::string PhyloTreeToNhx(const bpp::PhyloTree &tree);

  /// \brief Print a tree in an ostream in newick (parenthesis) format.
  /// \return Nothing but edit an output stream by writing the tree in parenthesis format (default as Newick).
  static void write(const bpp::PhyloTree& tree /// tree to print.
      , std::ostream& os /// stream object.
      , const Text_format output = Text_format::newick /// print tree in format.
      , const std::string& description = "" /// Tree description.
      , const bool reconciliation_mode = false /// replace bootstraps by node name.
  );

};


#endif //PHYLASOLVER_IO_H
