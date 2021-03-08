#pragma once

#include "BaseReconciliationModel.hpp"


class GTBaseReconciliationInterface: public BaseReconciliationModel {
  public:
  GTBaseReconciliationInterface(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &recModelInfo):
    BaseReconciliationModel(speciesTree,
        geneSpeciesMapping,
        recModelInfo)
  {}

  virtual ~GTBaseReconciliationInterface() {}
  virtual double computeLogLikelihood() = 0;
  virtual void setInitialGeneTree(PLLUnrootedTree &tree) = 0;
  virtual bool inferMLScenario(Scenario &scenario, 
      bool stochastic = false) = 0;
  virtual bool isParsimony() const = 0;
  virtual void setRoot(pll_unode_t * root) = 0;
  virtual pll_unode_t *getRoot() = 0;
  virtual void invalidateAllCLVs() = 0;
  virtual void invalidateCLV(unsigned int geneNodeIndex) = 0;
  virtual void enableMADRooting(bool enable)= 0;
  virtual pll_unode_t *computeMLRoot() = 0;


};


template <class REAL>
class GTBaseReconciliationModel: public GTBaseReconciliationInterface {
public:  
  GTBaseReconciliationModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const RecModelInfo &recModelInfo); 

  virtual ~GTBaseReconciliationModel() {}
  
  virtual double computeLogLikelihood();
  virtual void setInitialGeneTree(PLLUnrootedTree &tree);
  virtual bool inferMLScenario(Scenario &scenario, bool stochastic = false);
  virtual bool isParsimony() const {return false;}
  virtual void setRoot(pll_unode_t * root) {_geneRoot = root;}
  virtual pll_unode_t *getRoot() {return _geneRoot;}
  virtual void invalidateAllCLVs();
  virtual void invalidateCLV(unsigned int geneNodeIndex);

  virtual void updateCLV(pll_unode_t *geneNode) = 0;
  virtual void computeGeneRootLikelihood(pll_unode_t *virtualRoot) = 0;
  virtual REAL getGeneRootLikelihood(pll_unode_t *root) const = 0;
  virtual REAL getGeneRootLikelihood(pll_unode_t *root, pll_rnode_t *speciesRoot) = 0;
  // Called by inferMLScenario
  // fills scenario with the best likelihood set of events that 
  // would lead to the subtree of geneNode under speciesNode
  // Assumes that all the CLVs are filled
  virtual bool backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      Scenario &scenario,
      bool isVirtualRoot = false,
      bool stochastic = false);
  virtual void computeProbability(pll_unode_t *geneNode, 
      pll_rnode_t *speciesNode, 
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
  pll_unode_t *getLeft(pll_unode_t *node, bool virtualRoot) const;  
  pll_unode_t *getRight(pll_unode_t *node, bool virtualRoot) const;
  pll_unode_t *getGeneSon(pll_unode_t *node, bool left, bool virtualRoot = false) const;
  void updateCLVs(bool invalidate = true);
  virtual void enableMADRooting(bool enable);
  virtual pll_unode_t *computeMLRoot();

  virtual REAL getLikelihoodFactor() {return REAL(1.0);}

private:
  void mapGenesToSpecies();
  void computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot);
  void updateCLVsRec(pll_unode_t *node);
  void markInvalidatedNodes();
  void markInvalidatedNodesRec(pll_unode_t *node);
  virtual void computeLikelihoods();
  double getSumLikelihood();

protected:
  pll_unode_t *_geneRoot;
  std::vector<pll_rnode_t *> _geneToSpeciesLCA;
  // gene ids in postorder 
  std::vector<unsigned int> _geneIds;
  unsigned int _maxGeneId;
  // set of invalid CLVs. All the CLVs from these CLVs to
  // the root(s) need to be recomputed
  std::unordered_set<unsigned int> _invalidatedNodes;
  std::vector<bool> _isCLVUpdated;
  std::vector<pll_unode_t *> _allNodes;
  PLLUnrootedTree *_pllUnrootedTree;
  bool _madRootingEnabled;
  std::vector<double> _madProbabilities;
};


template <class REAL>
GTBaseReconciliationModel<REAL>::GTBaseReconciliationModel(
    PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const RecModelInfo &recModelInfo):
  GTBaseReconciliationInterface(speciesTree, 
      geneSpeciesMapping,
      recModelInfo),
  _geneRoot(0),
  _maxGeneId(1),
  _pllUnrootedTree(nullptr),
  _madRootingEnabled(false)
{
}


template <class REAL>
void GTBaseReconciliationModel<REAL>::initFromUtree(pll_utree_t *tree) {
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
void GTBaseReconciliationModel<REAL>::mapGenesToSpecies()
{
  this->_geneToSpecies.resize(_allNodes.size());
  for (auto node: _allNodes) {
    if (!node->next) {
      std::string speciesName = 
        this->_geneNameToSpeciesName[std::string(node->label)]; 
      this->_geneToSpecies[node->node_index] = 
        this->_speciesNameToId[speciesName];
    }
  }
  this->_numberOfCoveredSpecies = 0;
  this->_speciesCoverage = std::vector<unsigned int>(
      this->_allSpeciesNodesCount, 0);
  for (auto node: _allNodes) {
    if (!this->_speciesCoverage[this->_geneToSpecies[node->node_index]]) {
      this->_numberOfCoveredSpecies++;
    }
    this->_speciesCoverage[this->_geneToSpecies[node->node_index]]++;
  }
  this->onSpeciesTreeChange(nullptr);
}

template <class REAL>
void GTBaseReconciliationModel<REAL>::setInitialGeneTree(PLLUnrootedTree &tree)
{
  _pllUnrootedTree = &tree;
  initFromUtree(tree.getRawPtr());
  mapGenesToSpecies();
  _maxGeneId = static_cast<unsigned int>(_allNodes.size() - 1);
  _geneToSpeciesLCA.resize(_maxGeneId + 1);
  invalidateAllCLVs();
}

template <class REAL>
void GTBaseReconciliationModel<REAL>::getRoots(std::vector<pll_unode_t *> &roots,
    const std::vector<unsigned int> &geneIds)
{
  roots.clear();
  if (this->_info.rootedGeneTree && _geneRoot) {
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
double GTBaseReconciliationModel<REAL>::computeLogLikelihood()
{
  this->beforeComputeLogLikelihood();
  auto root = getRoot();
  updateCLVs();
  computeLikelihoods();
  if (this->_info.rootedGeneTree) {
    setRoot(computeMLRoot());
    while (root != getRoot()) {
      updateCLVs(false);
      computeLikelihoods();
      root = getRoot();
      setRoot(computeMLRoot());
    }
  }
  
  auto res = getSumLikelihood();
  return res;
}

template <class REAL>
pll_unode_t *GTBaseReconciliationModel<REAL>::getGeneSon(pll_unode_t *node, bool left, bool virtualRoot) const
{
  if (left) {
    return getLeft(node, virtualRoot);
  } else {
    return getRight(node, virtualRoot);
  }
}

template <class REAL>
pll_unode_t *GTBaseReconciliationModel<REAL>::getLeft(pll_unode_t *node, 
    bool virtualRoot) const
{
  return virtualRoot ? node->next : node->next->back;
}

template <class REAL>
pll_unode_t *GTBaseReconciliationModel<REAL>::getRight(pll_unode_t *node, 
    bool virtualRoot) const
{
  return virtualRoot ? node->next->back : node->next->next->back;
}

template <class REAL>
void GTBaseReconciliationModel<REAL>::markInvalidatedNodesRec(
    pll_unode_t *node)
{
  _isCLVUpdated[node->node_index] = false;
  if (node->back->next) {
    markInvalidatedNodesRec(node->back->next);
    markInvalidatedNodesRec(node->back->next->next);
  }
}

template <class REAL>
void GTBaseReconciliationModel<REAL>::markInvalidatedNodes()
{
  for (auto nodeIndex: _invalidatedNodes) {
    auto node = _allNodes[nodeIndex];
    markInvalidatedNodesRec(node);
  }
  _invalidatedNodes.clear();
}

template <class REAL>
void GTBaseReconciliationModel<REAL>::updateCLVsRec(pll_unode_t *node)
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
      _geneToSpeciesLCA[gid] = this->_speciesTree.getNode(
          this->_geneToSpecies[gid]);
    } else { // gene internal node
      auto left = currentNode->next->back->node_index;
      auto right = currentNode->next->next->back->node_index;
      _geneToSpeciesLCA[gid] = 
        this->_speciesTree.getLCA(this->_geneToSpeciesLCA[left], 
            this->_geneToSpeciesLCA[right]); 
    }

    updateCLV(currentNode);
    nodes.pop();
    _isCLVUpdated[currentNode->node_index] = true;
  }
}

template <class REAL>
void GTBaseReconciliationModel<REAL>::updateCLVs(bool invalidate)
{
  if (invalidate) {
    switch (this->_likelihoodMode) {
      case PartialLikelihoodMode::PartialGenes:
        this->invalidateAllSpeciesCLVs();
      break;
      case PartialLikelihoodMode::PartialSpecies:
        this->invalidateAllCLVs();
      break;
      case PartialLikelihoodMode::NoPartial:
        this->invalidateAllSpeciesCLVs();
        this->invalidateAllCLVs();
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
void GTBaseReconciliationModel<REAL>::invalidateCLV(unsigned int nodeIndex)
{
  _invalidatedNodes.insert(nodeIndex);
}

template <class REAL>
void GTBaseReconciliationModel<REAL>::invalidateAllCLVs()
{
  _isCLVUpdated = std::vector<bool>(_maxGeneId + 1, false);
}
 

template <class REAL>
void GTBaseReconciliationModel<REAL>::enableMADRooting(bool enable)
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
void GTBaseReconciliationModel<REAL>::computeMLRoot(pll_unode_t *&bestGeneRoot, pll_rnode_t *&bestSpeciesRoot) 
{
  std::vector<pll_unode_t *> roots;
  getRoots(roots, _geneIds);
  REAL max = isParsimony() ? 
    REAL(-std::numeric_limits<double>::infinity()) 
    : REAL();
  for (auto root: roots) {
    for (auto speciesNode: this->_allSpeciesNodes) {
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
pll_unode_t *GTBaseReconciliationModel<REAL>::computeMLRoot()
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
double GTBaseReconciliationModel<REAL>::getSumLikelihood()
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
    return log(total) - log(this->getLikelihoodFactor()); 
  } else {
    return double(total);
  }
}


template <class REAL>
void GTBaseReconciliationModel<REAL>::computeLikelihoods()
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

template <class REAL>
bool GTBaseReconciliationModel<REAL>::inferMLScenario(Scenario &scenario, bool stochastic)
{
  // make sure the CLVs are filled
  invalidateAllCLVs();
  updateCLVs();
  computeLikelihoods(); 
  auto ll = getSumLikelihood();
  assert(ll == 0.0 || (std::isnormal(ll) && ll <= 0.0));

  pll_unode_t *geneRoot = 0;
  pll_rnode_t *speciesRoot = 0;
  computeMLRoot(geneRoot, speciesRoot);
  assert(geneRoot);
  assert(speciesRoot);
  scenario.setGeneRoot(geneRoot);
  scenario.setSpeciesTree(this->_speciesTree.getRawPtr());
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = geneRoot->node_index + _maxGeneId + 1;
  scenario.setVirtualRootIndex(virtualRoot.node_index);
  scenario.initBlackList(_maxGeneId, this->_speciesTree.getNodesNumber());
  return  backtrace(&virtualRoot, speciesRoot, scenario, true, stochastic);
}

template <class REAL>
bool GTBaseReconciliationModel<REAL>::backtrace(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
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







