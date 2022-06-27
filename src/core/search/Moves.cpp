#include <search/Moves.hpp>
#include <trees/JointTree.hpp>
#include <IO/Logger.hpp>

// constants taken from RAXML
#define DEF_LH_EPSILON            0.1
#define OPT_LH_EPSILON            0.1
#define RAXML_PARAM_EPSILON       0.001  //0.01
#define RAXML_BFGS_FACTOR         1e7
#define RAXML_BRLEN_SMOOTHINGS    5
#define RAXML_BRLEN_DEFAULT       0.1
#define RAXML_BRLEN_MIN           1.0e-6
#define RAXML_BRLEN_MAX           100.
#define RAXML_BRLEN_TOLERANCE     1.0e-9
#define RAXML_FREERATE_MIN        0.001
#define RAXML_FREERATE_MAX        100.
#define RAXML_BRLEN_SCALER_MIN    0.01
#define RAXML_BRLEN_SCALER_MAX 100.
  


std::shared_ptr<SPRMove> SPRMove::createSPRMove(unsigned int pruneIndex, 
    unsigned int regraftIndex, 
    const std::vector<unsigned int> &path) {
  return std::make_shared<SPRMove>(pruneIndex, regraftIndex, path);
}


static void optimizeBranchesSlow(JointTree &tree,
    const std::vector<corax_unode_t *> &nodesToOptimize)
{
    auto root = tree.getTreeInfo()->root;
    // could be incremental and thus faster
    auto treeinfo = tree.getTreeInfo();
    auto ratecats = tree.getModel().num_ratecats();
    std::vector<unsigned int> params_indices(ratecats, 0);
    tree.computeLibpllLoglk(); // update CLVs
    for (unsigned int j = 0; j < 2; ++j) {
      for (unsigned int i = 0; i < nodesToOptimize.size(); ++i) {
          corax_treeinfo_set_root(treeinfo, nodesToOptimize[i]);
          double oldLoglk = tree.computeLibpllLoglk(true);
          double newLoglk = corax_opt_optimize_branch_lengths_local(
              treeinfo->partitions[0],
              treeinfo->root,
              &params_indices[0],
              RAXML_BRLEN_MIN,
              RAXML_BRLEN_MAX,
              RAXML_BRLEN_TOLERANCE,
              RAXML_BRLEN_SMOOTHINGS,
              0,
              true);
         assert(oldLoglk <= newLoglk);
      }
    }
    corax_treeinfo_set_root(treeinfo, root);
}


SPRMove::SPRMove(unsigned int pruneIndex, unsigned int regraftIndex, const std::vector<unsigned int> &path):
  _pruneIndex(pruneIndex),
  _regraftIndex(regraftIndex),
  _path(path)
{
}

std::shared_ptr<SPRRollback> SPRMove::applyMove(JointTree &tree)
{
  auto root = tree.getRoot();
  auto prune = tree.getNode(_pruneIndex);
  auto regraft = tree.getNode(_regraftIndex);
  assert (prune && prune->next);
  assert (regraft && regraft);
  tree.invalidateCLV(prune->next->back);
  tree.invalidateCLV(prune->next->next);
  tree.invalidateCLV(prune->next);
  tree.invalidateCLV(prune->next->next->back);
  tree.invalidateCLV(regraft);
  tree.invalidateCLV(regraft->back);
  for (auto branchIndex: _path) {
    tree.invalidateCLV(tree.getNode(branchIndex));
    tree.invalidateCLV(tree.getNode(branchIndex)->back);
  }
  corax_tree_rollback_t corax_rollback;
  std::vector<SavedBranch> savedBranches;
    
  _branchesToOptimize.push_back(prune);
  _branchesToOptimize.push_back(regraft->back);
  _branchesToOptimize.push_back(regraft);
  for (auto branchIndex: _path) {
    auto node = tree.getNode(branchIndex);
    _branchesToOptimize.push_back(node);
    if (_path.size() == 1) {
      _branchesToOptimize.push_back(node->next);
      _branchesToOptimize.push_back(node->next->next);
      node = node->back;
      _branchesToOptimize.push_back(node->next);
      _branchesToOptimize.push_back(node->next->next);
    }
  }
  for (auto branch: _branchesToOptimize) {
    savedBranches.push_back(branch);
  }
  assert(CORAX_SUCCESS == corax_utree_spr(prune, regraft, &corax_rollback));
  return std::make_shared<SPRRollback>(tree, corax_rollback, savedBranches, root);
}
  
void SPRMove::optimizeMove(JointTree &tree)
{
  optimizeBranchesSlow(tree, _branchesToOptimize);
  _branchesToOptimize.clear();
}

std::ostream& SPRMove::print(std::ostream & os) const {
  os << "SPR(";
  os << "prune:" <<_pruneIndex << ", regraft:" << _regraftIndex;
  os << ",_pathsize:" << _path.size();
  os << ")";
  return os;
}
  
void SPRMove::updatePath(JointTree &tree)
{
  auto prune = tree.getNode(_pruneIndex); 
  auto regraft = tree.getNode(_regraftIndex); 
  auto initialRegraft = regraft;
  auto initialPrune = prune;
  std::vector<corax_unode_t *> nodes;
  _path.clear();
  PLLUnrootedTree::orientTowardEachOther(&prune, &regraft, nodes);
  for (auto node: nodes) {
    if (node != regraft) {
      _path.push_back(node->node_index);
    }
  
  }
}

