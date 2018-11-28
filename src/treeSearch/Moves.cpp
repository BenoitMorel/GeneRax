#include <treeSearch/Moves.h>
#include <treeSearch/JointTree.h>
#include <Logger.hpp>
#include <Arguments.hpp>

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
  

std::shared_ptr<Move> Move::createNNIMove(int nodeIndex, bool left, bool blo, int bloRadius) {
  return make_shared<NNIMove>(nodeIndex, left, blo, bloRadius);
}

std::shared_ptr<Move> Move::createSPRMove(int pruneIndex, int regraftIndex, const vector<int> &path) {
  return make_shared<SPRMove>(pruneIndex, regraftIndex, path);
}

NNIMove::NNIMove(int nodeIndex, bool left, bool blo, int bloRadius):
    nodeIndex_(nodeIndex),
    left_(left),
    blo_(blo && !Arguments::noFelsensteinLikelihood),
    bloRadius_(bloRadius)
{
}

void NNIOptimizeBranches(JointTree &tree, pll_unode_t *edge, int radius)
{
    auto root = tree.getTreeInfo()->root;
    // could be incremental and thus faster
    unsigned int params_indices[4] = {0,0,0,0};
    auto treeinfo = tree.getTreeInfo();
    //for (int j = 0; j < 3; ++j) {
      pllmod_treeinfo_set_root(treeinfo.get(), edge);
      double oldLoglk = tree.computeLibpllLoglk(); // update CLVs
      double newLoglk = pllmod_opt_optimize_branch_lengths_local(
          treeinfo->partitions[0],
          edge,
          params_indices,
          RAXML_BRLEN_MIN,
          RAXML_BRLEN_MAX,
          RAXML_BRLEN_TOLERANCE,
          RAXML_BRLEN_SMOOTHINGS,
          radius,
          true);
   // }
    pllmod_treeinfo_set_root(treeinfo.get(), root);

}

void optimizeBranchesSlow(JointTree &tree,
    const std::vector<pll_unode_t *> &nodesToOptimize)
{
    auto root = tree.getTreeInfo()->root;
    // could be incremental and thus faster
    unsigned int params_indices[4] = {0,0,0,0};
    auto treeinfo = tree.getTreeInfo();
    for (int j = 0; j < 3; ++j) 
    for (unsigned int i = 0; i < nodesToOptimize.size(); ++i) {
        pllmod_treeinfo_set_root(treeinfo.get(), nodesToOptimize[i]);
        double oldLoglk = tree.computeLibpllLoglk(); // update CLVs
        double newLoglk = pllmod_opt_optimize_branch_lengths_local(
            treeinfo->partitions[0],
            treeinfo->root,
            params_indices,
            RAXML_BRLEN_MIN,
            RAXML_BRLEN_MAX,
            RAXML_BRLEN_TOLERANCE,
            RAXML_BRLEN_SMOOTHINGS,
            0,
            true);
       assert(oldLoglk <= newLoglk);
    }
    pllmod_treeinfo_set_root(treeinfo.get(), root);

}
void NNIMove::applyBLO(JointTree &tree, pll_unode_t *edge) {
  if (!Arguments::incr) {
    std::vector<pll_unode_t *> nodesToOptimize;
    nodesToOptimize.push_back(edge);
    if (edge->next) {
        nodesToOptimize.push_back(edge->next);
        nodesToOptimize.push_back(edge->next->next);
    }
    if (edge->back && edge->back->next) {
        nodesToOptimize.push_back(edge->back->next);
        nodesToOptimize.push_back(edge->back->next->next);
    }
    optimizeBranchesSlow(tree, nodesToOptimize);
  } else {
    NNIOptimizeBranches(tree, edge, bloRadius_);
  }
}

bool equals(pll_unode_t *node1, pll_unode_t *node2) {
  return node1 == node2 || (node1->next && (node1->next == node2 || node1->next->next == node2));
}

std::shared_ptr<Rollback> NNIMove::applyMove(JointTree &tree) {
    auto edge = tree.getNode(nodeIndex_);
    auto treeRoot = tree.getReconciliationEvaluation()->getRoot();
    bool rootAffected = false;
    if (treeRoot) {
      rootAffected = equals(treeRoot, edge) || equals(treeRoot, edge->back);
      if (edge->next) {
        rootAffected = rootAffected || equals(treeRoot, edge->next->back);
        rootAffected = rootAffected || equals(treeRoot, edge->next->next->back);
      }
      if (edge->back->next) {
        rootAffected = rootAffected || equals(treeRoot, edge->back->next->back);
        rootAffected = rootAffected || equals(treeRoot, edge->back->next->next->back);
      }
      if (rootAffected) {
        tree.getReconciliationEvaluation()->setRoot(0);
      }
    }
    unsigned int move = left_ ? PLL_UTREE_MOVE_NNI_LEFT : PLL_UTREE_MOVE_NNI_RIGHT;
    pll_tree_rollback_t rollback;
    assert(PLL_SUCCESS == pllmod_utree_nni(edge, move, &rollback));
    if (blo_) {
        applyBLO(tree, edge);
    }
    return std::make_shared<NNIRollback>(tree, rollback, rootAffected);
}

ostream& NNIMove::print(ostream & os) const {
  os << "NNI(" << nodeIndex_ << ",";
  os << (left_ ? "left":"right") << ",";
  os << ")";
  return os;
}


SPRMove::SPRMove(int pruneIndex, int regraftIndex, const vector<int> &path):
  pruneIndex_(pruneIndex),
  regraftIndex_(regraftIndex),
  path_(path)
{
}

void printNode(pll_unode_s *node)
{
  Logger::info << node;
  if (node->next) {
    Logger::info << " " << node->next << " " << node->next->next;
  }
  Logger::info << endl;
}



std::shared_ptr<Rollback> SPRMove::applyMove(JointTree &tree)
{
  bool blo = true;
  auto prune = tree.getNode(pruneIndex_);
  auto regraft = tree.getNode(regraftIndex_);
  pll_tree_rollback_t pll_rollback;

  vector<SavedBranch> savedBranches;;
    
  vector<pll_unode_t *> branchesToOptimize;
  if (blo) {
    branchesToOptimize.push_back(prune);
    branchesToOptimize.push_back(prune->next->back);
    branchesToOptimize.push_back(regraft->back);
    branchesToOptimize.push_back(regraft);
    for (int branchIndex: path_) {
      savedBranches.push_back(tree.getNode(branchIndex));
      branchesToOptimize.push_back(tree.getNode(branchIndex));
    }
    
  } 
  
  assert(PLL_SUCCESS == pllmod_utree_spr(prune, regraft, &pll_rollback));
  auto rollback = std::make_shared<SPRRollback>(tree, pll_rollback, savedBranches);
  optimizeBranchesSlow(tree, branchesToOptimize);
  return rollback;
}

ostream& SPRMove::print(ostream & os) const {
  os << "SPR(";
  os << "prune:" <<pruneIndex_ << ", regraft:" << regraftIndex_;
  os << ",path_size:" << path_.size();
  os << ")";
  return os;
}


