#include "DSTagger.hpp"
#include <IO/Logger.hpp>
#include <climits>
#include <sstream>


DSTagger::DSTagger(PLLUnrootedTree &tree):_tree(tree),
  _isRootDup(false),
  _clvs(_tree.getDirectedNodeNumber() + 3)
{
  for (auto node: _tree.getPostOrderNodes()) {
    auto &clv = _clvs[node->node_index];
    _tagNode(node, clv);
  }
  
  // here we create a fake node with fake children
  // to run the _tagNode function on virtual roots
  corax_unode_t fakeNode; // virtual root
  corax_unode_t fakeNext;
  corax_unode_t fakeNextNext;
  fakeNode.next = &fakeNext;
  fakeNode.next->next = &fakeNextNext;
  unsigned int bestScore = UINT_MAX;
  for (auto branch: _tree.getBranches()) {
    fakeNext.back = branch;
    fakeNextNext.back = branch->back;
    CLV clv;
    _tagNode(&fakeNode, clv);
    if (clv.score == bestScore) {
      _bestRoots.push_back(branch);
    } else if (clv.score < bestScore) {
      bestScore = clv.score;
      _bestRoots.clear();
      _bestRoots.push_back(branch);
      if (clv.isDup) {
        _isRootDup = true;
      }
    }
  }
  auto rootBranch = getRoot();
  _clvs[rootBranch->node_index].isRoot = true;
  _clvs[rootBranch->back->node_index].isRoot = true;
  

  for (unsigned int i = 0; i < 3; ++i) {
    _roots[i].node_index = _tree.getDirectedNodeNumber() + i;
    _roots[i].next = &_roots[(i + 1) % 3];
    _clvs[_roots[i].node_index].isDup = _isRootDup;
  }
  _roots[0].back = nullptr;
  _roots[1].back = rootBranch;
  _roots[2].back = rootBranch->back;

  _rootFromNode(&_roots[0]);
}
  
void DSTagger::fillUpTraversal(corax_unode_t *node,
      TaxaSet &set)
{
  auto parents = getSpeciationAncestorNodes(node);
  for (auto parent: parents) {
    auto parentUp = parent;
    orientUp(parentUp);
    assert(parentUp != parent);
    auto cousinsParent = getThirdNode(parent, parentUp);
    fillWithChildren(cousinsParent->back, set);  
  }
}

void DSTagger::_tagNode(corax_unode_t *node, CLV &clv)
{
  if (!node->next) {
    // leaf case, nothing to do
    clv.clade.insert(node->clv_index);
    return;
  }
  auto &leftCLV = _clvs[node->next->back->node_index];
  auto &rightCLV = _clvs[node->next->next->back->node_index];
  auto &leftClade = leftCLV.clade;
  auto &rightClade = rightCLV.clade;
  clv.score = leftCLV.score + rightCLV.score;
  clv.clade.insert(leftClade.begin(), leftClade.end());
  clv.clade.insert(rightClade.begin(), rightClade.end());
  clv.isDup = clv.clade.size() != 
    (leftClade.size() + rightClade.size());
  if (clv.isDup) {
    if (clv.clade == leftClade || clv.clade == rightClade) {
      if (leftClade == rightClade) {
        clv.score += 1;
      } else {
        clv.score += 2;
      }
    } else {
      clv.score += 3;
    }
  }
}
  
void DSTagger::_rootFromNode(corax_unode_t *node)
{
  auto &clv = _clvs[node->node_index];
  clv.up = node;
  
  if (node->next) {
    auto &clv2 = _clvs[node->next->node_index];
    auto &clv3 = _clvs[node->next->next->node_index];
    clv2.isRoot = clv.isRoot;
    clv3.isRoot = clv.isRoot;
    clv2.isDup = clv.isDup;
    clv3.isDup = clv.isDup;
    clv2.up = node;
    clv3.up = node;
    _rootFromNode(node->next->back);
    _rootFromNode(node->next->next->back);
  }
}


void DSTagger::print()
{
  auto root = getRoot();
  Logger::info << "DS Tagged:" << std::endl;
  TaggerUNodePrinter printer(_clvs);
  Logger::info << _tree.getNewickString(printer, 
      root,
      true) << std::endl;
}
    

void DSTagger::TaggerUNodePrinter::operator()(corax_unode_t *node, 
        std::stringstream &ss) const
{
  if (!node->next) {
    ss << node->label;
  } else {
    ss << (clvs[node->node_index].isDup ? "D" : "S"); 
  } 
  ss << ":" << node->length;
}
  

static void _fillWithChildren(corax_unode_t *node,
    TaxaSet &set)
{
  if (!node->next) {
    set.insert(node->clv_index);
  } else{
    _fillWithChildren(node->next->back, set);
    _fillWithChildren(node->next->next->back, set);
  }
}

void DSTagger::fillWithChildren(corax_unode_t *node,
      TaxaSet &set)
{
  orientUp(node);
  _fillWithChildren(node, set);  
}

std::vector<corax_unode_t *> 
DSTagger::getSpeciationAncestorNodes(corax_unode_t *node)
{
  std::vector<corax_unode_t *> ancestors;
  auto clv = &_clvs[node->node_index];
  while (true) { //!clv->isRoot) {
    orientUp(node);
    assert(node);
    if (clv->isRoot) {
      if (_isRootDup) {
        break;
      }
      if (_roots[1].back == node) {
        node = &_roots[1];
      } else if (_roots[2].back == node) {
        node = &_roots[2];
      } else {
        assert(false);
      }
    } else {
      node = node->back;
    }
    if (!node) { // we reached the root
      break;
    }
    clv = &_clvs[node->node_index];
    if (!clv->isDup) {
      ancestors.push_back(node);
    }
  }
  return ancestors;
}


static void _fillWithInternalDescendants(corax_unode_t *node,
    std::vector<corax_unode_t*> &descendants)
{
  if (node->next) {
    descendants.push_back(node);
    _fillWithInternalDescendants(node->next->back, descendants);  
    _fillWithInternalDescendants(node->next->next->back, descendants);
  }
}

void DSTagger::fillWithInternalDescendants(corax_unode_t *node,
      std::vector<corax_unode_t*> &descendants)
{
  if (node->next) {
    _fillWithInternalDescendants(node->next->back, descendants);  
    _fillWithInternalDescendants(node->next->next->back, descendants);
  }
}

