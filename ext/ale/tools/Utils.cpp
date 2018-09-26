// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 11/01/17.
//

#include "Utils.h"

#include <regex>

#include "PhyloTreeToolBox.h"

void Utils::printNodeContent(const bpp::PhyloTree& tree, const std::shared_ptr<bpp::PhyloNode>& node, std::ostream& os){
  ///Print Node Content (infos, father, sons).

  //First, print father and itself
  if(tree.getRoot() != node and tree.hasFather(node)) {
    auto father = tree.getFather(node);
    auto edgeFatherNode = tree.getEdgeLinking(father, node);
    if(edgeFatherNode) {
      if (edgeFatherNode->hasLength())
        os << father << " -(" << tree.getEdgeLinking(father, node)->getLength() << ")-> " << node << std::endl;
      else
        os << father << " ---> " << node << std::endl;
    }
    else
      os << father << " -x-> " << node << std::endl;

  } else {
    os << node << " (root)" << std::endl;
  }

  //Then print sons
  if(not tree.isLeaf(node)) {
    auto sons = tree.getSons(node);
    for (auto son: sons) {
      auto edgeNodeSon = tree.getEdgeLinking(node, son);
      if(edgeNodeSon) {
        if (edgeNodeSon->hasLength())
          os << "         '-(" << tree.getEdgeLinking(node, son)->getLength() << ")-> " << son << std::endl;
        else
          os << "         '---> " << son << std::endl;
      }
      else
        os << "         '-x-> " << son << std::endl;
    }
  }
}

void Utils::print_temp(bpp::PhyloTree &tree, std::ostream& os) {
  /// Print all bifurcations in a given tree.
  PhyloTreeToolBox::applyInPOT(tree,
                     [&tree, &os](std::shared_ptr<bpp::PhyloNode>& node)
                      {
                        printNodeContent(tree, node, os);
                      });
}

bool Utils::stringIsNumber(const std::string &s) {
  /// Check if a string is an int
  if(s.empty()) return false;
  char* p;
  std::strtol(s.c_str(), &p, 10);
  return (*p == 0);
}

bool Utils::string_comp(const std::string &a, const std::string &b, const bool case_sensitive) {
  if (a.length()==b.length()) {
    return std::equal(b.begin(), b.end(),
                      a.begin(),
                      case_sensitive ? Utils::case_sensitive_char_comp : Utils::case_insensitive_char_comp);
  }
  else {
    return false;
  }
}

std::vector<std::string> Utils::splitString(const std::string &str, const char *c, const bool concatenate_tokens) {
  /// Split a string according to delimiters defined in c.
  std::string buff{""};
  std::vector<std::string> ll;

  for (auto n:str) {
    // find if n is in the delimiters list in const char* c
    bool n_is_delimiter = false;
    std::size_t i = 0;

    while (!n_is_delimiter && c[i] != '\0') {
      n_is_delimiter = (n == c[i++]);
    }

    if (!n_is_delimiter){
      buff += n;
    }else if (n_is_delimiter && buff != "") {
      ll.push_back(buff);
      buff = "";
    } else if(n_is_delimiter && not concatenate_tokens) {
      ll.push_back(buff);
    }
  }
  if (buff != "") ll.push_back(buff);

  return ll;
}

bool Utils::isYes(const char *optarg) {
  /// Check if a given const char* means "yes" or "true".
  return
      string_comp(optarg, "true", false)
      or string_comp(optarg, "yes", false)
      or string_comp(optarg, "t", false)
      or string_comp(optarg, "y", false);
}

bool Utils::strmatch_regex(const std::string &str, const std::string &pattern) {
  auto regex = std::regex(pattern);
  return std::regex_match(str, regex);
}

unsigned int Utils::count(const std::string &str, const std::string &sub) {
  if (sub.length() == 0 or sub.length() > str.length()) return 0;
  unsigned int count = 0;
  for (std::size_t offset = str.find(sub); offset != std::string::npos;
       offset = str.find(sub, offset + sub.length()))
  {
    ++count;
  }
  return count;
}

std::string Utils::extractFilename(const std::string &str) {
  // First, extract the correct file name (without the path). So, find the last path delimiter (which is the beginning of the filename)
  // and keep only the string after.
  std::size_t last_path_delimiter = str.find_last_of(PATH_SLASH);

  std::size_t output_name_begin = 0;
  if(last_path_delimiter != std::string::npos) output_name_begin = last_path_delimiter + 1;

  // Remove the input file name extension.
  std::size_t last_dot = str.find_last_of('.');
  if(last_dot == std::string::npos or output_name_begin > last_dot) {
    last_dot = str.size();
  }

  if(output_name_begin < last_dot)
    return str.substr(output_name_begin, last_dot - output_name_begin);

  return str;
}

std::string Utils::replace(const std::string& str
                           , const std::string& old_substr
                           , const std::string& new_substr
) {
  return std::regex_replace(str, std::regex(old_substr), new_substr);
}

std::string Utils::replace(const std::string& str
                           , const char old_char
                           , const char new_char
) {
  std::string res(str);
  for(char& c: res) {
    if(c == old_char)
      c = new_char;
  }
  return res;
}

std::ostream &operator<<(std::ostream &os, const bpp::PhyloTree &tree) {
  /// Print a bpp::PhyloTree in newick format.
  bpp::Newick newick = bpp::Newick();
  newick.write(tree, os);
  return os;
}

std::ostream &operator<<(std::ostream &os, const std::shared_ptr<bpp::PhyloNode> &node_ptr) {
  /// Print the name of a std::shared_ptr<bpp::PhyloNode>.
  if(node_ptr) {
    if (node_ptr->hasName())
      os << node_ptr->getName();
    else
      os << "Noname";
  } else
    os << "nullptr";
  return os;
}
