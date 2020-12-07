#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/GeneSpeciesMapping.hpp>
#include <util/Scenario.hpp>
#include <IO/Logger.hpp>
#include <util/enums.hpp>
#include <cmath>
#include <unordered_set>
#include <maths/ScaledValue.hpp>
#include <trees/PLLRootedTree.hpp>
#include <maths/Random.hpp>



//#define IS_PROBA(x) ((REAL(0.0) <= REAL(x) && REAL(x)  <= REAL(1.0)))
//#define ASSERT_PROBA(x) assert(IS_PROBA(x));
#define IS_PROBA(x) true
#define ASSERT_PROBA(x) assert(true);

double log(ScaledValue v);


typedef std::vector< std::vector <double> > RatesVector;

/**
 *  Interface and common implementations for 
 *  all the reconciliation likelihood computation
 *  classes
 */
class ReconciliationModelInterface {
public:
  virtual ~ReconciliationModelInterface() {}
  
  virtual void setInitialGeneTree(PLLUnrootedTree &tree) = 0;
  
  /*
   * Set the per-species lineage rates
   */
  virtual void setRates(const RatesVector &rates) = 0;
  
  /**
   * (incrementally) compute and return the likelihood of the gene tree 
   */
  virtual double computeLogLikelihood() = 0;

  /**
   * If implemented, rollback to the state before the
   * last "computeLogLikelihood(false)" call
   */
  virtual void rollbackToLastState() = 0;

  virtual bool isParsimony() const = 0;

  /**
   *  Get/set the root of the gene tree (only relevant in rooted gene tree mode)
   */ 
  virtual void setRoot(pll_unode_t * root) = 0;
  virtual pll_unode_t *getRoot() = 0;
  
  /**
   * Compute and set the maximum likelihood root. Relevant in both rooted and unrooted gene tree modes
   */
  virtual pll_unode_t *computeMLRoot() = 0;

  virtual void enableMADRooting(bool enable) = 0;

  virtual void setPartialLikelihoodMode(PartialLikelihoodMode mode) = 0;
  
  /**
   * CLV invalidation for partial likelihood computation
   */
  virtual void invalidateAllCLVs() = 0;
  virtual void invalidateCLV(unsigned int geneNodeIndex) = 0;
  virtual void invalidateAllSpeciesCLVs() = 0;
  /**
   *  Fill scenario with the maximum likelihood set of 
   *  events that would lead to the  current tree
   **/
  virtual bool inferMLScenario(Scenario &scenario, bool stochastic = false) = 0;

  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate) = 0;

  virtual void setFractionMissingGenes(const std::string &fractionMissingFile) = 0;

  /**
   * If return true, the evaluated likelihood is the product of
   * the likelihoods of all possible originations. Else, the only
   * possible origination is at the species tree root.
   */
  virtual bool sumOverAllOriginations() const = 0;
};


/*
 * Introduce the REAL template (normal double or infinite precision) 
 * Implement util methods shared by its children
 */
template <class REAL>
class AbstractReconciliationModel: public ReconciliationModelInterface {
public:
  AbstractReconciliationModel(const AbstractReconciliationModel &) = delete;
  AbstractReconciliationModel & operator = (const AbstractReconciliationModel &) = delete;
  AbstractReconciliationModel(AbstractReconciliationModel &&) = delete;
  AbstractReconciliationModel & operator = (AbstractReconciliationModel &&) = delete;
  
  
  // overload from parent
  AbstractReconciliationModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &recModelInfo); 
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  virtual ~AbstractReconciliationModel() {}
  // overload from parent
  virtual void setRates(const RatesVector &rates) = 0;
  // overload from parent
  virtual double computeLogLikelihood();
  // overload from parent 
  virtual bool isParsimony() const {return false;}
  // overload from parent 
  virtual void rollbackToLastState() {assert(false);}
  // overload from parent
  virtual void setRoot(pll_unode_t * root) {_geneRoot = root;}
  // overload from parent
  virtual pll_unode_t *getRoot() {return _geneRoot;}
  // overload from parent
  virtual void invalidateAllCLVs();
  // overload from parent
  virtual void invalidateCLV(unsigned int geneNodeIndex);
  // overload from parent
  virtual void invalidateAllSpeciesCLVs() {_allSpeciesNodesInvalid = true;}
  // overload from parent
  virtual bool inferMLScenario(Scenario &scenario, bool stochastic = false);
  // overload from parent
  virtual void setPartialLikelihoodMode(PartialLikelihoodMode mode) {_likelihoodMode = mode;};
  // overload from parents
  virtual void setFractionMissingGenes(const std::string &fractionMissingFile);
  // overload from parents
  virtual bool sumOverAllOriginations() const {return true;}
protected:
  // called by the constructor
  virtual void initSpeciesTree();
  // Called by computeLogLikelihood
  virtual void updateCLV(pll_unode_t *geneNode) = 0;
  // Called by computeLogLikelihood
  virtual void computeGeneRootLikelihood(pll_unode_t *virtualRoot) = 0;
  // Called by computeLogLikelihood
  virtual REAL getGeneRootLikelihood(pll_unode_t *root) const = 0;
  virtual REAL getGeneRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) = 0;
  virtual REAL getLikelihoodFactor() const = 0;
  virtual void recomputeSpeciesProbabilities() = 0;
  // Called by inferMLScenario
  // fills scenario with the best likelihood set of events that 
  // would lead to the subtree of geneNode under speciesNode
  // Can assume that all the CLVs are filled
  virtual bool backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false,
      bool stochastic = false);
 
  virtual void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      REAL &proba,
      bool isVirtualRoot = false,
      Scenario *scenario = nullptr,
      Scenario::Event *event = nullptr,
      bool stochastic = false) = 0;
  
  void initFromUtree(pll_utree_t *tree);
  /**
   *  - In rooted gene tree mode, and if the gene tree already has a virtual root,
   *  return this root and its direct neighbors
   *  - Else, return all the possible virtual roots
   */
  void getRoots(std::vector<pll_unode_t *> &roots,
    const std::vector<unsigned int> &geneIds);

  /**
   *  Get the left or right child of node. If node is a virtual 
   *  root, the implementation is different
   */
  pll_unode_t *getLeft(pll_unode_t *node, bool virtualRoot) const;  
  pll_unode_t *getRight(pll_unode_t *node, bool virtualRoot) const;
  pll_unode_t *getGeneSon(pll_unode_t *node, bool left, bool virtualRoot = false) const;
  pll_unode_t *getLeftRepeats(pll_unode_t *node, bool virtualRoot);  
  pll_unode_t *getRightRepeats(pll_unode_t *node, bool virtualRoot); 

  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);
  
  void updateCLVs(bool invalidate = true);
  virtual void enableMADRooting(bool enable);
  virtual pll_unode_t *computeMLRoot();
protected:
  RecModelInfo _info;
  pll_unode_t *_geneRoot;
  bool _mlEventProba;
  unsigned int _allSpeciesNodesCount;
  std::vector <pll_rnode_t *> _speciesNodesToUpdate;
  std::vector <pll_rnode_t *> _allSpeciesNodes;
  std::vector<unsigned int> _geneToSpecies;
  // gene ids in postorder 
  std::vector<unsigned int> _geneIds;
  std::vector<pll_rnode_t *> _geneToSpeciesLCA;
  unsigned int _maxGeneId;
  PartialLikelihoodMode _likelihoodMode;
  std::vector<double> _fm;
  PLLRootedTree &_speciesTree;
  virtual void beforeComputeLogLikelihood(); 
  virtual void afterComputeLogLikelihood() {};
  pll_rnode_t *getSpeciesSon(pll_rnode_t *node, bool left) {return left ? getSpeciesLeft(node) : getSpeciesRight(node);}
  pll_rnode_t *getSpeciesLeft(pll_rnode_t *node) {return _speciesLeft[node->node_index];}
  pll_rnode_t *getSpeciesRight(pll_rnode_t *node) {return _speciesRight[node->node_index];}
  pll_rnode_t *getSpeciesParent(pll_rnode_t *node) {return _speciesParent[node->node_index];}
  pll_rnode_t *getPrunedRoot() {return _prunedRoot;}
private:
  void mapGenesToSpecies();
  void computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot);
  virtual void computeLikelihoods();
  double getSumLikelihood();
  void updateCLVsRec(pll_unode_t *node);
  void markInvalidatedNodes();
  void markInvalidatedNodesRec(pll_unode_t *node);
  bool fillPrunedNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd = nullptr);  
  
  
  std::map<std::string, std::string> _geneNameToSpeciesName;
  std::map<std::string, unsigned int> _speciesNameToId;
  std::vector<unsigned int> _speciesCoverage;
  unsigned int _numberOfCoveredSpecies;
  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  std::unordered_set<unsigned int> _invalidatedNodes;
  
  std::unordered_set<pll_rnode_t *> _invalidatedSpeciesNodes;
  bool _allSpeciesNodesInvalid;

  // is the CLV up to date?
  std::vector<bool> _isCLVUpdated;
  std::vector<pll_unode_t *> _allNodes;
 
  // left, right and parent species vectors, 
  // index with the species nodex_index
  std::vector<pll_rnode_t *> _speciesLeft;
  std::vector<pll_rnode_t *> _speciesRight;
  std::vector<pll_rnode_t *> _speciesParent;
  pll_rnode_t *_prunedRoot;
  PLLUnrootedTree *_pllUnrootedTree;
  bool _madRootingEnabled;
  std::vector<double> _madProbabilities;

};


  
template <class REAL>
AbstractReconciliationModel<REAL>::AbstractReconciliationModel(
    PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &recModelInfo):
  _info(recModelInfo),
  _geneRoot(0),
  _mlEventProba(false),
  _maxGeneId(1),
  _likelihoodMode(PartialLikelihoodMode::PartialGenes),
  _speciesTree(speciesTree),
  _geneNameToSpeciesName(geneSpeciesMapping.getMap()),
  _allSpeciesNodesInvalid(true),
  _pllUnrootedTree(nullptr),
  _madRootingEnabled(false)
{
  initSpeciesTree();
  setFractionMissingGenes(_info.fractionMissingFile);
}

template <class REAL>
void AbstractReconciliationModel<REAL>::initFromUtree(pll_utree_t *tree) {
  auto treeSize = tree->tip_count + tree->inner_count;
  auto nodesNumber = tree->tip_count + 3 * tree->inner_count;
  _geneIds.clear();
  _allNodes.resize(nodesNumber);
  std::vector<bool> marked(nodesNumber, false);
  for (unsigned int i = 0; i < treeSize; ++i) {
    auto node = tree->nodes[i];
    _allNodes[node->node_index] = node;
    _geneIds.push_back(node->node_index);
    if  (node->next) {
      node = node->next;
      _allNodes[node->node_index] = node;
      _geneIds.push_back(node->node_index);
      node = node->next;
      _allNodes[node->node_index] = node;
      _geneIds.push_back(node->node_index);
    }
  }

}


template <class REAL>
void AbstractReconciliationModel<REAL>::mapGenesToSpecies()
{
  _geneToSpecies.resize(_allNodes.size());
  for (auto node: _allNodes) {
    if (!node->next) {
      std::string speciesName = _geneNameToSpeciesName[std::string(node->label)]; 
      _geneToSpecies[node->node_index] = _speciesNameToId[speciesName];
    }
  }
  _numberOfCoveredSpecies = 0;
  _speciesCoverage = std::vector<unsigned int>(_allSpeciesNodesCount, 0);
  for (auto node: _allNodes) {
    if (!_speciesCoverage[_geneToSpecies[node->node_index]]) {
      _numberOfCoveredSpecies++;
    }
    _speciesCoverage[_geneToSpecies[node->node_index]]++;
  }

  onSpeciesTreeChange(nullptr);
}

template <class REAL>
void AbstractReconciliationModel<REAL>::setInitialGeneTree(PLLUnrootedTree &tree)
{
  _pllUnrootedTree = &tree;
  initFromUtree(tree.getRawPtr());
  mapGenesToSpecies();
  _maxGeneId = static_cast<unsigned int>(_allNodes.size() - 1);
  _geneToSpeciesLCA.resize(_maxGeneId + 1);
  invalidateAllCLVs();
}
  
template <class REAL>
bool AbstractReconciliationModel<REAL>::fillPrunedNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd)
{
  bool addMyself = true;
  if (nodesToAdd) {
    addMyself = (nodesToAdd->find(node) != nodesToAdd->end());
  }
  if (getSpeciesLeft(node)) {
    assert(getSpeciesRight(node));
    addMyself |= fillPrunedNodesPostOrder(getSpeciesLeft(node), nodes, nodesToAdd);
    addMyself |= fillPrunedNodesPostOrder(getSpeciesRight(node), nodes, nodesToAdd);
  }
  if (addMyself) {
    nodes.push_back(node);
  }
  return addMyself;
}


static bool fillNodesPostOrder(pll_rnode_t *node, 
    std::vector<pll_rnode_t *> &nodes, 
    std::unordered_set<pll_rnode_t *> *nodesToAdd = nullptr)  
{
  bool addMyself = true;
  if (nodesToAdd) {
    addMyself = (nodesToAdd->find(node) != nodesToAdd->end());
  }
  if (node->left) {
    assert(node->right);
    addMyself |= fillNodesPostOrder(node->left, nodes, nodesToAdd);
    addMyself |= fillNodesPostOrder(node->right, nodes, nodesToAdd);
  }
  if (addMyself) {
    nodes.push_back(node);
  }
  return addMyself;
}


template <class REAL>
void AbstractReconciliationModel<REAL>::initSpeciesTree()
{
  _allSpeciesNodesCount = _speciesTree.getNodesNumber();
  _speciesLeft = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesRight = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesParent = std::vector<pll_rnode_t *>(_allSpeciesNodesCount, nullptr);
  _speciesNameToId.clear();
  onSpeciesTreeChange(nullptr);
  for (auto node: _allSpeciesNodes) {
    if (!node->left) {
      _speciesNameToId[node->label] = node->node_index;
    }
  }
}

template <class REAL>
void AbstractReconciliationModel<REAL>::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  if (!nodesToInvalidate) {
    _allSpeciesNodesInvalid = true;
  } else {
    assert(nodesToInvalidate->size());
    for (auto node: *nodesToInvalidate) {
      while (node) {
        _invalidatedSpeciesNodes.insert(node);
        node = getSpeciesParent(node);
      }
    }
    _invalidatedSpeciesNodes.insert(nodesToInvalidate->begin(), nodesToInvalidate->end());
  }
  _allSpeciesNodes.clear();
  fillNodesPostOrder(_speciesTree.getRoot(), _allSpeciesNodes);
  for (auto speciesNode: _allSpeciesNodes) {
    auto e = speciesNode->node_index;
    _speciesLeft[e] = speciesNode->left;
    _speciesRight[e] = speciesNode->right;
    _speciesParent[e] = speciesNode->parent;
  }
  _prunedRoot = _speciesTree.getRoot();
  if (_info.pruneSpeciesTree && _speciesCoverage.size()) {
    std::vector<pll_rnode_t *> pruned(_allSpeciesNodesCount, nullptr);
    for (auto speciesNode: _allSpeciesNodes) {
      auto e = speciesNode->node_index;
      if (!speciesNode->left) {
        pruned[e] = (_speciesCoverage[e] ? speciesNode : nullptr);   
      } else {
        auto left = _speciesLeft[e];
        auto right = _speciesRight[e];
        auto prunedLeft = pruned[left->node_index];
        auto prunedRight = pruned[right->node_index];
        if (prunedLeft && prunedRight) {
          _speciesLeft[e] = prunedLeft;
          _speciesRight[e] = prunedRight;
          pruned[e] = speciesNode;
          _speciesParent[prunedLeft->node_index] = speciesNode;
          _speciesParent[prunedRight->node_index] = speciesNode;
          _prunedRoot = speciesNode;
        } else if (!prunedLeft && prunedRight) {
          pruned[e] = prunedRight;
        } else if (prunedLeft && !prunedRight) {
          pruned[e] = prunedLeft;
        }
      }
    }
    _allSpeciesNodes.clear();
    fillPrunedNodesPostOrder(getPrunedRoot(), _allSpeciesNodes);
  }
  assert(_allSpeciesNodes.size()); // && _allSpeciesNodes.back() == _speciesTree.getRoot());
}


template <class REAL>
void AbstractReconciliationModel<REAL>::beforeComputeLogLikelihood()
{
  if (_allSpeciesNodesInvalid) { // update everything
    _speciesNodesToUpdate = _allSpeciesNodes;
  } else if (_invalidatedSpeciesNodes.size()) { // partial update
    // here, fill _speciesNodesToUpdate with the invalid nodes
    _speciesNodesToUpdate.clear();
    fillPrunedNodesPostOrder(getPrunedRoot(), _speciesNodesToUpdate, &_invalidatedSpeciesNodes);
  } else {
    _speciesNodesToUpdate.clear();
  } 
  _allSpeciesNodesInvalid = false;
  _invalidatedSpeciesNodes.clear();
  //assert(!_speciesNodesToUpdate.size() || _speciesNodesToUpdate.back() == getPrunedRoot());
  recomputeSpeciesProbabilities();
}

template <class REAL>
void AbstractReconciliationModel<REAL>::getRoots(std::vector<pll_unode_t *> &roots,
    const std::vector<unsigned int> &geneIds)
{
  roots.clear();
  if (_info.rootedGeneTree && _geneRoot) {
    roots.push_back(_geneRoot);
    if (_geneRoot->next) {
      roots.push_back(_geneRoot->next);
      roots.push_back(_geneRoot->next->next);
    }
    if (_geneRoot->back->next) {
      roots.push_back(_geneRoot->back->next);
      roots.push_back(_geneRoot->back->next->next);
    }
    return;
  }
  std::vector<bool> marked(geneIds.size(), false);
  for (auto id: geneIds) {
    auto node = _allNodes[id];
    if (marked[node->node_index] || marked[node->back->node_index]) {
      continue;
    }
    roots.push_back(node->back);
    marked[node->node_index] = true;
  }
}
  
template <class REAL>
double AbstractReconciliationModel<REAL>::computeLogLikelihood()
{
  beforeComputeLogLikelihood();
  auto root = getRoot();
  updateCLVs();
  computeLikelihoods();
  if (_info.rootedGeneTree) {
    setRoot(computeMLRoot());
    while (root != getRoot()) {
      updateCLVs(false);
      computeLikelihoods();
      root = getRoot();
      setRoot(computeMLRoot());
    }
  }
  
  auto res = getSumLikelihood();
  afterComputeLogLikelihood();
  return res;
}

template <class REAL>
pll_unode_t *AbstractReconciliationModel<REAL>::getGeneSon(pll_unode_t *node, bool left, bool virtualRoot) const
{
  if (left) {
    return getLeft(node, virtualRoot);
  } else {
    return getRight(node, virtualRoot);
  }
}

template <class REAL>
pll_unode_t *AbstractReconciliationModel<REAL>::getLeft(pll_unode_t *node, bool virtualRoot) const
{
  return virtualRoot ? node->next : node->next->back;
}

template <class REAL>
pll_unode_t *AbstractReconciliationModel<REAL>::getRight(pll_unode_t *node, bool virtualRoot) const
{
  return virtualRoot ? node->next->back : node->next->next->back;
}

template <class REAL>
pll_unode_t *AbstractReconciliationModel<REAL>::getLeftRepeats(pll_unode_t *node, bool virtualRoot)
{
  return getLeft(node, virtualRoot);
}

template <class REAL>
pll_unode_t *AbstractReconciliationModel<REAL>::getRightRepeats(pll_unode_t *node, bool virtualRoot)
{
  return getRight(node, virtualRoot);
}

template <class REAL>
void AbstractReconciliationModel<REAL>::markInvalidatedNodesRec(pll_unode_t *node)
{
  _isCLVUpdated[node->node_index] = false;
  if (node->back->next) {
    markInvalidatedNodesRec(node->back->next);
    markInvalidatedNodesRec(node->back->next->next);
  }
}

template <class REAL>
void AbstractReconciliationModel<REAL>::markInvalidatedNodes()
{
  for (auto nodeIndex: _invalidatedNodes) {
    auto node = _allNodes[nodeIndex];
    markInvalidatedNodesRec(node);
  }
  _invalidatedNodes.clear();
}

template <class REAL>
void AbstractReconciliationModel<REAL>::updateCLVsRec(pll_unode_t *node)
{
  if (_isCLVUpdated[node->node_index]) {
    return;
  }
  std::stack<pll_unode_t *> nodes;
  nodes.push(node);
  while (!nodes.empty()) {
    auto currentNode = nodes.top();
    if (currentNode->next) {
      bool waitForChildren = false;
      auto left = getLeft(currentNode, false);
      auto right = getRight(currentNode, false);
      if (!_isCLVUpdated[left->node_index]) {
        nodes.push(left);
        waitForChildren = true;
      }
      if (!_isCLVUpdated[right->node_index]) {
        nodes.push(right);
        waitForChildren = true;
      }
      if (waitForChildren) {
        continue;
      }
    }
    
    // update LCA
    auto gid = currentNode->node_index;
    if (!currentNode->next) { // gene leaf
      _geneToSpeciesLCA[gid] = _speciesTree.getNode(_geneToSpecies[gid]);
    } else { // gene internal node
      auto left = currentNode->next->back->node_index;
      auto right = currentNode->next->next->back->node_index;
      _geneToSpeciesLCA[gid] = 
        this->_speciesTree.getLCA(_geneToSpeciesLCA[left], _geneToSpeciesLCA[right]); 
    }

    updateCLV(currentNode);
    nodes.pop();
    _isCLVUpdated[currentNode->node_index] = true;
  }
}

template <class REAL>
void AbstractReconciliationModel<REAL>::updateCLVs(bool invalidate)
{
  if (invalidate) {
    switch (_likelihoodMode) {
      case PartialLikelihoodMode::PartialGenes:
        invalidateAllSpeciesCLVs();
      break;
      case PartialLikelihoodMode::PartialSpecies:
        invalidateAllCLVs();
      break;
      case PartialLikelihoodMode::NoPartial:
        invalidateAllSpeciesCLVs();
        invalidateAllCLVs();
      break;
    }
    markInvalidatedNodes();
  }
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  for (auto root: roots) {
    updateCLVsRec(root);
    updateCLVsRec(root->back);
  }
}

template <class REAL>
void AbstractReconciliationModel<REAL>::invalidateCLV(unsigned int nodeIndex)
{
  _invalidatedNodes.insert(nodeIndex);
}
  
template <class REAL>
void AbstractReconciliationModel<REAL>::setFractionMissingGenes(const std::string &fractionMissingFile)
{
  _fm = std::vector<double>(_allSpeciesNodesCount, 0.0);
  if (!fractionMissingFile.size()) {
    return;
  }
  std::ifstream is(fractionMissingFile);
  std::string species;
  double fm;
  while (is >> species >> fm) {
    _fm[_speciesNameToId[species]] = fm;
  }
}
  
template <class REAL>
void AbstractReconciliationModel<REAL>::invalidateAllCLVs()
{
  _isCLVUpdated = std::vector<bool>(_maxGeneId + 1, false);
}
 

template <class REAL>
void AbstractReconciliationModel<REAL>::enableMADRooting(bool enable)
{
  _madRootingEnabled = enable;
  double beta = -2;
  double epsilon = 0.01;
  if (enable) {
    _madProbabilities = _pllUnrootedTree->getMADRelativeDeviations();
    double sum = 0.0;
    for (auto &v: _madProbabilities) {
      v += epsilon;
      v = pow(v, beta);
      sum += v;
    }
    for (auto &v: _madProbabilities) {
      v /= sum;
    }
  } 
}

template <class REAL>
void AbstractReconciliationModel<REAL>::computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot) 
{
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  REAL max = isParsimony() ? 
    REAL(-std::numeric_limits<double>::infinity()) 
    : REAL();
  assert(sumOverAllOriginations());
  for (auto root: roots) {
    for (auto speciesNode: _allSpeciesNodes) {
      REAL ll = getGeneRootLikelihood(root, speciesNode);
      if (_madRootingEnabled) {
        ll *= _madProbabilities[root->node_index];
      }
      if (max < ll) {
        max = ll;
        bestGeneRoot = root;
        bestSpeciesRoot = speciesNode;
      }
    }
  }
}

template <class REAL>
pll_unode_t *AbstractReconciliationModel<REAL>::computeMLRoot()
{
  pll_unode_t *bestRoot = 0;
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  REAL max = isParsimony() ? 
    REAL(-std::numeric_limits<double>::infinity()) 
    : REAL();
  for (auto root: roots) {
    REAL rootProba = getGeneRootLikelihood(root);
    if (_madRootingEnabled) {
      rootProba *= _madProbabilities[root->node_index];
    }
    if (max < rootProba) {
      bestRoot = root;
      max = rootProba;
    }
  }
  assert(bestRoot);
  return bestRoot;
}

template <class REAL>
double AbstractReconciliationModel<REAL>::getSumLikelihood()
{
  REAL total = REAL();
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  //Logger::info << "HEY" << std::endl;
  if (!isParsimony()) {
    for (auto root: roots) {
      auto ll = getGeneRootLikelihood(root);
      if (_madRootingEnabled) {
        ll *= _madProbabilities[root->node_index];
        //Logger::info << root->node_index << " " << _madProbabilities[root->node_index] << std::endl;
      }

      total += ll;
    }
  } else {
    total =  REAL(-std::numeric_limits<double>::infinity()); 
    for (auto root: roots) {
      auto v = getGeneRootLikelihood(root);
      if (total < v) {
        total = v;
      }
    }
  }
  bool applyLog = !isParsimony();
  if (applyLog) {
    return log(total) - log(getLikelihoodFactor()); 
  } else {
    return double(total);
  }
}


template <class REAL>
void AbstractReconciliationModel<REAL>::computeLikelihoods()
{
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  for (auto root: roots) {
    pll_unode_t virtualRoot;
    virtualRoot.next = root;
    virtualRoot.node_index = root->node_index + _maxGeneId + 1;
    computeGeneRootLikelihood(&virtualRoot);
  }
}

 
static std::string speciesTreeString(pll_rnode_t *node)
{
  if (!node->left) {
    return std::to_string(node->node_index);
  } 
  std::string res = "(";
  res += speciesTreeString(node->left);
  res += ",";
  res += speciesTreeString(node->right);
  res += ")";
  res += std::to_string(node->node_index);
  return res;
}

static std::string geneTreeString(pll_unode_t *node, const std::vector<unsigned int> toSpecies)
{
  std::string res;
  if (!node->next) {
    //res += std::to_string(node->node_index);
    //res += "-";
    res += std::to_string(toSpecies[node->node_index]);
  } else {
    auto left = node->next->back;
    auto right = node->next->next->back;
    res += "(";
    res += geneTreeString(left, toSpecies);
    res += ",";
    res += geneTreeString(right, toSpecies);
    res += ")";
    //res += std::to_string(node->node_index);
  }
  return res;
}

template <class REAL>
bool AbstractReconciliationModel<REAL>::inferMLScenario(Scenario &scenario, bool stochastic)
{
  // make sure the CLVs are filled
  invalidateAllCLVs();
  updateCLVs();
  computeLikelihoods(); 
  auto ll = getSumLikelihood();
  assert(ll == 0.0 || (std::isnormal(ll) && ll <= 0.0));

  pll_unode_t *geneRoot = 0;
  pll_rnode_t *speciesRoot = 0;
  auto ll1 = computeLogLikelihood();
  if(!stochastic) {
    auto saveMLEventProba = _mlEventProba;
    _mlEventProba = true;
    auto ll2 = computeLogLikelihood();
    _mlEventProba = saveMLEventProba;
  }
  computeMLRoot(geneRoot, speciesRoot);
  assert(geneRoot);
  assert(speciesRoot);
  scenario.setGeneRoot(geneRoot);
  scenario.setSpeciesTree(_speciesTree.getRawPtr());
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = geneRoot->node_index + _maxGeneId + 1;
  scenario.setVirtualRootIndex(virtualRoot.node_index);
  scenario.initBlackList(_maxGeneId, _speciesTree.getNodesNumber());
  return  backtrace(&virtualRoot, speciesRoot, scenario, true, stochastic);
}
  

static pll_unode_t *getOther(pll_unode_t *ref, pll_unode_t *n1, pll_unode_t *n2)
{
  return (ref == n1) ? n2 : n1;
}

template <class REAL>
bool AbstractReconciliationModel<REAL>::backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot, 
      bool stochastic) 
{
  pll_unode_t *leftGeneNode = 0;     
  pll_unode_t *rightGeneNode = 0;     
  bool isGeneLeaf = !geneNode->next;
  if (!isGeneLeaf) {
    leftGeneNode = this->getLeft(geneNode, isVirtualRoot);
    rightGeneNode = this->getRight(geneNode, isVirtualRoot);
  }
  REAL temp;
  Scenario::Event event;
  computeProbability(geneNode, speciesNode, temp, isVirtualRoot, &scenario, &event, stochastic);
  scenario.addEvent(event);
  bool ok = true;
  // safety check
  switch(event.type) {
  case ReconciliationEventType::EVENT_S:
    if (!event.cross) {
      ok &= backtrace(leftGeneNode, this->getSpeciesLeft(speciesNode), scenario, false, stochastic); 
      ok &= backtrace(rightGeneNode, this->getSpeciesRight(speciesNode), scenario, false, stochastic); 
    } else {
      ok &= backtrace(leftGeneNode, this->getSpeciesRight(speciesNode), scenario, false, stochastic); 
      ok &= backtrace(rightGeneNode, this->getSpeciesLeft(speciesNode), scenario, false, stochastic); 
    }
    break;
  case ReconciliationEventType::EVENT_D:
    ok &= backtrace(leftGeneNode, speciesNode, scenario, false, stochastic); 
    ok &= backtrace(rightGeneNode, speciesNode, scenario, false, stochastic); 
    break;
  case ReconciliationEventType::EVENT_SL:
    ok &= backtrace(geneNode, event.pllDestSpeciesNode, scenario, isVirtualRoot, stochastic); 
    break;
  case ReconciliationEventType::EVENT_T:
    ok &= backtrace(event.pllTransferedGeneNode, event.pllDestSpeciesNode, scenario, false, stochastic);
    ok &= backtrace(getOther(event.pllTransferedGeneNode, leftGeneNode, rightGeneNode),
        speciesNode, scenario, false, stochastic);
    break;
  case ReconciliationEventType::EVENT_TL:
    ok &= backtrace(geneNode, event.pllDestSpeciesNode, scenario, isVirtualRoot, stochastic);
    break;
  case ReconciliationEventType::EVENT_None:
    break;
  default:
    ok  = false;
    break;
  }
  return ok;
}


template<class REAL> 
REAL getRandom(REAL max)
{
  return max * Random::getProba();
}
template<class C, class REAL>
int sampleIndex(const C &container)
{
  REAL sum = REAL();
  for (auto value: container) {
    sum += value;
  }
  if (sum == REAL()) {
    return -1;
  }
  REAL stopAt = getRandom(sum);
  sum = REAL();
  unsigned int index = 0;
  for (auto value: container) {
    sum += value;
    if (stopAt < sum) {
      return index;
    }
    index++;
  }
  assert(false);
}

