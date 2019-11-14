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



//#define IS_PROBA(x) ((x) >= REAL(0) && (x) <= REAL(1))
//#define ASSERT_PROBA(x) assert(IS_PROBA(x));
#define ASSERT_PROBA(x) assert(true);




/**
 *  Interface and common implementations for 
 *  all the reconciliation likelihood computation
 *  classes
 */
class ReconciliationModelInterface {
public:
  virtual ~ReconciliationModelInterface() {}
  
  virtual void setInitialGeneTree(pll_utree_t *tree) = 0;
  
  /*
   * Set the per-species lineage rates
   *  @param dupRate
   *  @param lossRate
   *  @param transferRate
   */
  virtual void setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates) = 0;
  
  /**
   * (incrementally) compute and return the likelihood of the gene tree 
   */
  virtual double computeLogLikelihood() = 0;
  
  /**
   *  Get/set the root of the gene tree (only relevant in rooted gene tree mode)
   */ 
  virtual void setRoot(pll_unode_t * root) = 0;
  virtual pll_unode_t *getRoot() = 0;
  
  /**
   * Compute and set the maximum likelihood root. Relevant in both rooted and unrooted gene tree modes
   */
  virtual pll_unode_t *computeMLRoot() = 0;
  

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
  virtual void inferMLScenario(Scenario &scenario) = 0;

  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate) = 0;
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
      bool rootedGeneTree);
  virtual void setInitialGeneTree(pll_utree_t *tree);
  virtual ~AbstractReconciliationModel() {}
  // overload from parent
  virtual void setRates(const std::vector<double> &dupRates,
      const std::vector<double> &lossRates,
      const std::vector<double> &transferRates) = 0;
  // overload from parent
  virtual double computeLogLikelihood();
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
  virtual void inferMLScenario(Scenario &scenario);
  // overload from parent
  virtual void setPartialLikelihoodMode(PartialLikelihoodMode mode) {_likelihoodMode = mode;};
protected:
  // called by the constructor
  virtual void initSpeciesTree();
  // Called when computeLogLikelihood is called for the first time
  // Called by computeLogLikelihood
  virtual void updateCLV(pll_unode_t *geneNode) = 0;
  // Called by computeLogLikelihood
  virtual void computeRootLikelihood(pll_unode_t *virtualRoot) = 0;
  // Called by computeLogLikelihood
  virtual REAL getRootLikelihood(pll_unode_t *root) const = 0;
  virtual REAL getRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) = 0;
  virtual REAL getLikelihoodFactor() const = 0;
  // Called by inferMLScenario
  // fills scenario with the best likelihood set of events that 
  // would lead to the subtree of geneNode under speciesNode
  // Can assume that all the CLVs are filled
  virtual void backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false) = 0;
 
  
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
  pll_unode_t *getLeftRepeats(pll_unode_t *node, bool virtualRoot);  
  pll_unode_t *getRightRepeats(pll_unode_t *node, bool virtualRoot); 

  virtual void onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate);
  
  void updateCLVs();
  virtual pll_unode_t *computeMLRoot();
protected:
  pll_unode_t *_geneRoot;
  unsigned int _allSpeciesNodesCount;
  std::vector <pll_rnode_t *> _speciesNodesToUpdate;
  std::vector <pll_rnode_t *> _allSpeciesNodes;
  std::vector<unsigned int> _geneToSpecies;
  // gene ids in postorder 
  std::vector<unsigned int> _geneIds;
  unsigned int _maxGeneId;
private:
  void mapGenesToSpecies();
  void computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot);
  virtual void computeLikelihoods();
  double getSumLikelihood();
  void updateCLVsRec(pll_unode_t *node);
  void markInvalidatedNodes();
  void markInvalidatedNodesRec(pll_unode_t *node);
  void beforeComputeLogLikelihood(); 
  
  
  bool _rootedGeneTree;
  PLLRootedTree &_speciesTree;
  std::map<std::string, std::string> geneNameToSpeciesName_;
  std::map<std::string, unsigned int> speciesNameToId_;
  
  PartialLikelihoodMode _likelihoodMode;
  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  std::unordered_set<unsigned int> _invalidatedNodes;
  
  std::unordered_set<pll_rnode_t *> _invalidatedSpeciesNodes;
  bool _allSpeciesNodesInvalid;

  // is the CLV up to date?
  std::vector<bool> _isCLVUpdated;
  std::vector<pll_unode_t *> _allNodes;
  
};


  
template <class REAL>
AbstractReconciliationModel<REAL>::AbstractReconciliationModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    bool rootedGeneTree):
  _geneRoot(0),
  _maxGeneId(1),
  _rootedGeneTree(rootedGeneTree),
  _speciesTree(speciesTree),
  geneNameToSpeciesName_(geneSpeciesMapping.getMap()),
  _likelihoodMode(PartialLikelihoodMode::PartialGenes),
  _allSpeciesNodesInvalid(true)
{
  initSpeciesTree();
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
      std::string speciesName = geneNameToSpeciesName_[std::string(node->label)]; 
      _geneToSpecies[node->node_index] = speciesNameToId_[speciesName];
    }
  }
}

template <class REAL>
void AbstractReconciliationModel<REAL>::setInitialGeneTree(pll_utree_t *tree)
{
  initFromUtree(tree);
  mapGenesToSpecies();
  _maxGeneId = static_cast<unsigned int>(_allNodes.size() - 1);
  invalidateAllCLVs();
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
  beforeComputeLogLikelihood(); // todobenoit rename this function
  speciesNameToId_.clear();
  for (auto node: _allSpeciesNodes) {
    if (!node->left) {
      speciesNameToId_[node->label] = node->node_index;
    }
  }
}

template <class REAL>
void AbstractReconciliationModel<REAL>::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  if (!nodesToInvalidate) {
    _allSpeciesNodesInvalid = true;
  } else {
    _invalidatedSpeciesNodes.insert(nodesToInvalidate->begin(), nodesToInvalidate->end());
  }
}


template <class REAL>
void AbstractReconciliationModel<REAL>::beforeComputeLogLikelihood()
{
  if (_allSpeciesNodesInvalid) { // update everything
    //Logger::info << "All nodes invalid " << std::endl;
    _allSpeciesNodes.clear();
    fillNodesPostOrder(_speciesTree.getRoot(), _allSpeciesNodes);
    _speciesNodesToUpdate = _allSpeciesNodes;
  } else if (_invalidatedSpeciesNodes.size()) { // partial update
    // here, fill _speciesNodesToUpdate with the invalid nodes
    //Logger::info << "Some nodes invalid " << std::endl;
    _allSpeciesNodes.clear();
    fillNodesPostOrder(_speciesTree.getRoot(), _allSpeciesNodes);
    fillNodesPostOrder(_speciesTree.getRoot(), _speciesNodesToUpdate, &_invalidatedSpeciesNodes);
    _speciesNodesToUpdate = _allSpeciesNodes; // TO REMOVE todobenoit
  } // else: nothing to update
  _allSpeciesNodesInvalid = false;
  _invalidatedSpeciesNodes.clear();
  assert(_allSpeciesNodes.size() && _allSpeciesNodes.back() == _speciesTree.getRoot());
  assert(!_speciesNodesToUpdate.size() || _speciesNodesToUpdate.back() == _speciesTree.getRoot());
  //Logger::info << "all species: " << _allSpeciesNodes.size() << std::endl;
  //Logger::info << "species to update: " << _speciesNodesToUpdate.size() << std::endl;
}

template <class REAL>
void AbstractReconciliationModel<REAL>::getRoots(std::vector<pll_unode_t *> &roots,
    const std::vector<unsigned int> &geneIds)
{
  roots.clear();
  if (_rootedGeneTree && _geneRoot) {
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
  if (_rootedGeneTree) {
    setRoot(computeMLRoot());
    while (root != getRoot()) {
      updateCLVs();
      computeLikelihoods();
      root = getRoot();
      setRoot(computeMLRoot());
    }
  }
  return getSumLikelihood();
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
    updateCLV(currentNode);
    nodes.pop();
    _isCLVUpdated[currentNode->node_index] = true;
  }
}

template <class REAL>
void AbstractReconciliationModel<REAL>::updateCLVs()
{
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
void AbstractReconciliationModel<REAL>::invalidateAllCLVs()
{
  _isCLVUpdated = std::vector<bool>(_maxGeneId + 1, false);
}

template <class REAL>
void AbstractReconciliationModel<REAL>::computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot) 
{
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  REAL max = REAL();
  for (auto root: roots) {
    for (auto speciesNode: _allSpeciesNodes) {
      REAL ll = getRootLikelihood(root, speciesNode);
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
  REAL max = REAL();
  for (auto root: roots) {
    REAL rootProba = getRootLikelihood(root);
    if (max < rootProba) {
      bestRoot = root;
      max = rootProba;
    }
  }
  return bestRoot;
}

template <class REAL>
double AbstractReconciliationModel<REAL>::getSumLikelihood()
{
  REAL total = REAL();
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  for (auto root: roots) {
    total += getRootLikelihood(root);
  }
  return log(total) - log(getLikelihoodFactor()); 
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
    computeRootLikelihood(&virtualRoot);
  }
}
  
template <class REAL>
void AbstractReconciliationModel<REAL>::inferMLScenario(Scenario &scenario)
{
  // make sure the CLVs are filled
  invalidateAllCLVs();
  updateCLVs();
  computeLikelihoods(); 
  auto ll = getSumLikelihood();
  assert(std::isnormal(ll) && ll < 0.0);

  pll_unode_t *geneRoot = 0;
  pll_rnode_t *speciesRoot = 0;
  computeMLRoot(geneRoot, speciesRoot);
  assert(geneRoot);
  assert(speciesRoot);
  scenario.setGeneRoot(geneRoot);
  scenario.setSpeciesTree(_speciesTree.getRawPtr());
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = geneRoot->node_index + _maxGeneId + 1;
  scenario.setVirtualRootIndex(virtualRoot.node_index);
  backtrace(&virtualRoot, speciesRoot, scenario, true);

}
  

