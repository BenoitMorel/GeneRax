#include "LibpllParsers.hpp"
#include <fstream>
#include <streambuf>

pll_utree_t *LibpllParsers::readNewickFromFile(const string &newickFilename)
{
  ifstream t(newickFilename);
  if (!t)
    throw LibpllException("Could not load open newick file ", newickFilename);
  string str((istreambuf_iterator<char>(t)),
                       istreambuf_iterator<char>());
  return readNewickFromStr(str);
}


pll_utree_t *LibpllParsers::readNewickFromStr(const string &newickString)
{
  auto utree =  pll_utree_parse_newick_string_unroot(newickString.c_str());
  if (!utree) 
    throw LibpllException("Error while reading tree from string: ", newickString);
  return utree;
}

pll_rtree_t *LibpllParsers::readRootedFromFile(const string &newickFile)
{
  auto tree = pll_rtree_parse_newick(newickFile.c_str());
  if (!tree) {
    throw LibpllException("Error while reading tree fropm file: ", newickFile);
  }
  return tree;
}

