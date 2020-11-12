#include "PolytomySolver.hpp"
#include <trees/PLLRootedTree.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <limits>

const double BL_THRESHOLD = 0.000001;


  
struct Entry {
  unsigned int m1;
  unsigned int m2;
  unsigned int ym;
  unsigned int nbs;
  unsigned int dups;
  unsigned int losses;
  unsigned int ac; // ancestral copies
  std::vector<std::string> subtrees;
  Entry():
    m1(std::numeric_limits<unsigned int>::max()),
    m2(std::numeric_limits<unsigned int>::max()),
    ym(std::numeric_limits<unsigned int>::max()),
    nbs(0),
    dups(0),
    losses(0),
    ac(0)
  {}

  unsigned int getMinCost(unsigned int k)
  {
    if (k < m1) {
      return ym  + m1 - k;
    } else if (k > m2) {
      return ym + k - m2;
    } else {
      return ym;
    }
  }
};


static void computeCup(unsigned int l1,
    unsigned int l2,
    unsigned int yl,
    unsigned int r1,
    unsigned int r2,
    unsigned int yr,
    unsigned int &m1,
    unsigned int &m2,
    unsigned int &ym)
{
  m1 = 99999;
  if (l1 < r1) {
    if (l2 < r1) {
      // case 1
      ym = yl + yr + r1 - l2;
      m1 = l2;
      m2 = r1;
      return;
    } else if (r1 <= l2 && l2 <= r2) {
      // case 2
      ym = yl + yr;
      m1 = r1;
      m2 = l2;
      return;
    } else if (l2 > r2) {
      // case 3
      ym = yl + yr;
      m1 = r1;
      m2 = r2;
      return;
    } 
  } else if (r1 <= l1 && l1 <= r2) {
    if (l1 <= r2) {
      // case 4
      ym = yl + yr;
      m1 = l1;
      m2 = l2;
      return;
    } else if (l2 > r2) {
      // case 5
      ym = yl + yr;
      m1 = l1;
      m2 = r2;
      return;
    }
  } else if (l1 > r2 && l2 > r2) {
    ym = yl + yr + l1 - r2;
    m1 = r2;
    m2 = l1;
    return;
  }
}

void computeDupLoss(pll_rnode_t *node,
    std::vector<Entry> &entries,
    unsigned int k)
{
  auto &c = entries[node->node_index];
  if (!node->left) {
    // LEAF
    if (k >= c.nbs) {
      c.losses = k - c.nbs;
    } else {
      c.dups = c.nbs -k;
    }
  } else {
    // INTERNAL NODE
    auto cl = entries[node->left->node_index];
    auto cr = entries[node->right->node_index];
    
    if (k > c.nbs && 
        c.getMinCost(k) == 
        (cl.getMinCost(k-c.nbs) + cr.getMinCost(k-c.nbs))) {
      // case 1
      computeDupLoss(node->left, entries, k - c.nbs); 
      computeDupLoss(node->right, entries, k - c.nbs); 
    } else if (k < c.m1) {
      c.dups = c.m1 -k;
      computeDupLoss(node->left, entries, c.m1 - c.nbs); 
      computeDupLoss(node->right, entries, c.m1 - c.nbs); 
    } else if (k > c.m2) {
      c.losses = k - c.m2;
      computeDupLoss(node->left, entries, c.m2 - c.nbs); 
      computeDupLoss(node->right, entries, c.m2 - c.nbs); 
    }
  }
}

static void computeAncestralCopies(pll_rnode_t *node,
    std::vector<Entry> &entries,
    unsigned int parentsCopies)
{
  auto &c = entries[node->node_index];
  c.ac = parentsCopies + c.dups - c.losses;
  if (node->left) {
    computeAncestralCopies(node->left, entries, c.ac - c.nbs);
    computeAncestralCopies(node->right, entries, c.ac - c.nbs);
  }
}

std::string joinNodes(const std::string &n1, const std::string &n2) 
{
  if (!n1.size()) {
    return n2;
  } else if (!n2.size()) {
      return n1;   
  } else {
    return "(" + n1 + "," + n2 + ")";
  }
}

void PolytomySolver::solveSimpleInterface(
      PLLRootedTree &speciesTree,
      std::map<std::string, unsigned int> &speciesLabelsToSolve
      )
{
  unsigned int speciesNumber = speciesTree.getNodesNumber(); 
  std::vector<Entry> cells(speciesNumber);
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto spid = node->node_index;
    auto &c = cells[spid];
    std::string label = node->label;
    auto &m1 = c.m1;
    auto &m2 = c.m2;
    auto &ym = c.ym;
    
    if (speciesLabelsToSolve.find(label) != speciesLabelsToSolve.end()) {
      c.nbs = speciesLabelsToSolve[label];
    }
    if (!node->left) {
      // LEAF CASE
      if (c.nbs == 0) {
        m1 = m2 = 0;
        ym = 0;
      } else {
        m1 = m2 = c.nbs;
        ym = 0;
      }
    } else {
      auto &cl = cells[node->left->node_index];
      auto &cr = cells[node->right->node_index];
      computeCup(
          cl.m1, cl.m2, cl.ym,
          cr.m1, cr.m2, cr.ym,
          m1, m2, ym);
      m1 += c.nbs;
      m2 += c.nbs;
      //assert (m1 > c.nbs); 
      //assert (m2 > c.nbs); 
    }
    //std::cout << label << "\t" << m1 << "\t" << m2 << "\t" << ym << std::endl;
  }
  computeDupLoss(speciesTree.getRoot(), cells, 1);
  computeAncestralCopies(speciesTree.getRoot(), cells, 1);
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto spid = node->node_index;
    auto &c = cells[spid];
    //std::cout << node->label << "\tdup=" << c.dups << "\tloss=" << c.losses << "\tnbs=" << c.nbs << "\tac=" << c.ac << std::endl;
  }
  
  for (auto node: speciesTree.getPostOrderNodes()) {
    auto spid = node->node_index;
    auto &c = cells[spid];
    if (!node->left) {
      // LEAF
      std::string label(node->label);
      c.subtrees = std::vector<std::string>(c.nbs - c.dups, label);
      for (unsigned int i = 0; i < c.dups; ++i) {
        auto &t = c.subtrees[i%c.subtrees.size()];
        t = joinNodes(t, label);
      }
    } else {
      const auto &lTrees = cells[node->left->node_index].subtrees;
      const auto &rTrees = cells[node->right->node_index].subtrees;
      auto lSize = lTrees.size();
      auto rSize = rTrees.size();
      std::vector<std::string> specTrees(std::max(lSize, rSize));
      unsigned int speciations = std::min(lSize, rSize);
      // add speciations
      for (unsigned int i = 0; i < speciations; ++i) {
        specTrees[i] = joinNodes(lTrees[i], rTrees[i]);
      }
      // add left spec-loss
      for (unsigned int i = speciations; i < lSize; ++i) {
        specTrees[i] = lTrees[i];
      }
      // add right spec-loss
      for (unsigned int i = speciations; i < rSize; ++i) {
        specTrees[i] = rTrees[i];
      }
      // add ancestral nodes corresponding to current species
      for (unsigned int i = 0; i < c.nbs; ++i) {
        specTrees.push_back(std::string(node->label));
      }
      // now let's build the subtree with duplications
      c.subtrees.push_back("");
      for (unsigned int i = 0; i < c.dups + 1; ++i) {
        auto s = joinNodes(c.subtrees[0], specTrees[i]);
        c.subtrees[0] = s;
      }
      for (unsigned int i = c.dups + 1; i < specTrees.size(); ++i) {
        c.subtrees.push_back(specTrees[i]);
      }
    }
  }
  auto root = speciesTree.getRoot();
  auto &rootCell = cells[root->node_index];
  std::cout << root->label << " ";
  for (auto &s: rootCell.subtrees) {
    std::cout << s << " ";
  }
  std::cout << std::endl;



}

