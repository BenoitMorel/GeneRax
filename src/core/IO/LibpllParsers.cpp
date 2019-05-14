#include "LibpllParsers.hpp"
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <ParallelContext.hpp>

pll_utree_t *LibpllParsers::readNewickFromFile(const std::string &newickFilename)
{
  std::ifstream t(newickFilename);
  if (!t)
    throw LibpllException("Could not load open newick file ", newickFilename);
  
  std::string str((istreambuf_iterator<char>(t)),
                       istreambuf_iterator<char>());
  
  pll_utree_t *res = 0;
  try {
    res = readNewickFromStr(str);
  } catch (...) {
    throw LibpllException("Error while reading tree from file: ", newickFilename);
  }
  return res;
}


pll_utree_t *LibpllParsers::readNewickFromStr(const std::string &newickString)
{
  auto utree =  pll_utree_parse_newick_string_unroot(newickString.c_str());
  if (!utree) 
    throw LibpllException("Error while reading tree from std::string: ", newickString);
  return utree;
}

pll_rtree_t *LibpllParsers::readRootedFromFile(const std::string &newickFile)
{
  auto tree = pll_rtree_parse_newick(newickFile.c_str());
  if (!tree) {
    throw LibpllException("Error while reading tree fropm file: ", newickFile);
  }
  return tree;
}
  
void LibpllParsers::saveUtree(pll_unode_t *utree, 
  const std::string &fileName, 
  bool append)
{
  std::ofstream os(fileName, (append ? std::ofstream::app : std::ofstream::out));
  char *newick = pll_utree_export_newick_rooted(utree, 0);
  os << newick;
  os.close();
  free(newick);
}

std::vector<int> LibpllParsers::parallelGetTreeSizes(const std::vector<FamiliesFileParser::FamilyInfo> &families) 
{
  int treesNumber = families.size();
  std::vector<int> localTreeSizes((treesNumber - 1 ) / ParallelContext::getSize() + 1, 0);
  for (int i = ParallelContext::getBegin(treesNumber); i < ParallelContext::getEnd(treesNumber); i ++) {
    pll_utree_t *tree = LibpllParsers::readNewickFromFile(families[i].startingGeneTree);
    int taxa = tree->tip_count;
    localTreeSizes[i - ParallelContext::getBegin(treesNumber)] = taxa;
    pll_utree_destroy(tree, 0);
  }
  std::vector<int> treeSizes;
  ParallelContext::concatenateIntVectors(localTreeSizes, treeSizes);
  treeSizes.erase(remove(treeSizes.begin(), treeSizes.end(), 0), treeSizes.end());
  assert(treeSizes.size() == families.size());
  return treeSizes;
}

