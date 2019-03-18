#pragma once

extern "C" {
#include <pll.h>
}
#include <vector>
#include <string>
#include <IO/FamiliesFileParser.hpp>
using namespace std;

class LibpllException: public exception {
public:
  LibpllException(const string &s): msg_(s) {}
  LibpllException(const string &s1, 
      const string s2): msg_(s1 + s2) {}
  virtual const char* what() const throw() { return msg_.c_str(); }
  void append(const string &str) {msg_ += str;}

private:
  string msg_;
};

class LibpllParsers {
public:
  static pll_utree_t *readNewickFromFile(const string &newickFile);
  static pll_utree_t *readNewickFromStr(const string &newickSTring);
  static pll_rtree_t *readRootedFromFile(const string &newickFile);
  static vector<int> parallelGetTreeSizes(const vector<FamiliesFileParser::FamilyInfo> &families);
  static void saveUtree(pll_unode_t *utree, 
    const string &fileName, 
    bool append = false);
};

