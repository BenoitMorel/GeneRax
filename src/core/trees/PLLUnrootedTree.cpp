#include "PLLUnrootedTree.hpp"

#include <IO/LibpllParsers.hpp>
#include <IO/Logger.hpp>


void utreeDestroy(pll_utree_t *utree) {
  if(!utree)
    return;
  free(utree->nodes);
  free(utree);
}

static pll_utree_t *buildUtree(const std::string &str, bool isFile)
{
  if (isFile) {
    return LibpllParsers::readNewickFromFile(str);
  } else {
    return LibpllParsers::readNewickFromStr(str);
  }
}

PLLUnrootedTree::PLLUnrootedTree(const std::string &str, bool isFile):
  _tree(buildUtree(str, isFile), utreeDestroy)
{
  std::cerr << "Todo: check the nodes deallocation" << std::endl;
}

PLLUnrootedTree::PLLUnrootedTree(const std::vector<const char*> &labels,
    unsigned int seed):
  _tree(pllmod_utree_create_random(static_cast<unsigned int>(labels.size()), &labels[0], seed), utreeDestroy)
{

}

void PLLUnrootedTree::save(const std::string &fileName)
{
  LibpllParsers::saveUtree(_tree->nodes[0], fileName, false);
}

void PLLUnrootedTree::setMissingBranchLengths(double minBL)
{
  for (auto node: getLeaves()) {
    std::cerr << node->label << std::endl;
  }
  for (unsigned int i = 0; i < _tree->tip_count; ++i)
    if (0.0 == _tree->nodes[i]->length)
      _tree->nodes[i]->length = minBL;
  for (unsigned int i = _tree->tip_count; i < _tree->tip_count + _tree->inner_count; ++i) {
    if (0.0 == _tree->nodes[i]->length)
      _tree->nodes[i]->length = minBL;
    if (0.0 == _tree->nodes[i]->next->length)
      _tree->nodes[i]->next->length = minBL;
    if (0.0 == _tree->nodes[i]->next->next->length)
      _tree->nodes[i]->next->next->length = minBL;
  }  
}
  
CArrayRange<pll_unode_t*> PLLUnrootedTree::getNodes()
{
  return CArrayRange<pll_unode_t*>(_tree->nodes, getNodesNumber());
}

CArrayRange<pll_unode_t*> PLLUnrootedTree::getLeaves()
{
  return CArrayRange<pll_unode_t*>(_tree->nodes, getLeavesNumber());
}

CArrayRange<pll_unode_t*> PLLUnrootedTree::getInnerNodes()
{
  return CArrayRange<pll_unode_t*>(_tree->nodes + getLeavesNumber(), getNodesNumber());
}

unsigned int PLLUnrootedTree::getNodesNumber() const
{
  return getLeavesNumber() + getInnerNodesNumber();
}

unsigned int PLLUnrootedTree::getLeavesNumber() const
{
  return _tree->tip_count;
}

unsigned int PLLUnrootedTree::getInnerNodesNumber() const
{
  return _tree->inner_count * 3;
}
  
pll_unode_t *PLLUnrootedTree::getNode(unsigned int node_index)
{
  return _tree->nodes[node_index];
}

pll_unode_t *PLLUnrootedTree::getAnyInnerNode()
{
  return getNode(getNodesNumber() - 1);
}

