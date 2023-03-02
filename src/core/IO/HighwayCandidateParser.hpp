#pragma once
#include <string>
#include <trees/PLLRootedTree.hpp>

struct Highway {
  Highway():src(nullptr), dest(nullptr), proba(0.0) {}

  Highway(corax_rnode_t *src, corax_rnode_t *dest):
    src(src),
    dest(dest),
    proba(0.1)
  {}
  corax_rnode_t *src;
  corax_rnode_t *dest;
  double proba;
  
  friend std::ostream& operator<<(std::ostream& os, const Highway &h) {
    os << "(" << h.src->label << " -> " << h.dest->label << ")";
    return os;
  }
};
  

/**
 *  Functions to parse a file containing a list of candidate 
 *  transfer hightways to test. Example file:
 *
 *  from1,to1
 *  from2,to2
 *  from3,to3
 */
class HighwayCandidateParser {
public:
  static std::vector<Highway> parse(const std::string &candidateFile,
      PLLRootedTree &speciesTree);


};
