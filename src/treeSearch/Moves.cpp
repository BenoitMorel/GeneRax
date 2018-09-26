#include <treeSearch/Moves.h>
#include <treeSearch/JointTree.h>

// constants taken from RAXML
#define DEF_LH_EPSILON            0.1
#define OPT_LH_EPSILON            0.1
#define RAXML_PARAM_EPSILON       0.001  //0.01
#define RAXML_BFGS_FACTOR         1e7
#define RAXML_BRLEN_SMOOTHINGS    32
#define RAXML_BRLEN_DEFAULT       0.1
#define RAXML_BRLEN_MIN           1.0e-6
#define RAXML_BRLEN_MAX           100.
#define RAXML_BRLEN_TOLERANCE     1.0e-7
#define RAXML_FREERATE_MIN        0.001
#define RAXML_FREERATE_MAX        100.
#define RAXML_BRLEN_SCALER_MIN    0.01
#define RAXML_BRLEN_SCALER_MAX 100.
  
void Rollback::applyRollback() {
  assert(PLL_SUCCESS == pllmod_tree_rollback(&rollback_));
  tree_.updateBPPTree();
}


std::shared_ptr<Move> Move::createNNIMove(int nodeIndex, bool left, bool blo) {
  return make_shared<NNIMove>(nodeIndex, left, blo);
}

std::shared_ptr<Move> Move::createSPRMove(int pruneIndex, int regraftIndex) {
  return make_shared<SPRMove>(pruneIndex, regraftIndex);
}

NNIMove::NNIMove(int nodeIndex, bool left, bool blo):
    nodeIndex_(nodeIndex),
    left_(left),
    blo_(blo)
{
}

void optimizeBranches(JointTree &tree,
    const std::vector<pll_unode_t *> &nodesToOptimize)
{
    // could be incremental and thus faster
    unsigned int params_indices[4] = {0,0,0,0};
    auto treeinfo = tree.getTreeInfo();
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

}

void NNIMove::applyBLO(JointTree &tree, pll_unode_t *edge) {
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
    optimizeBranches(tree, nodesToOptimize);
}

std::shared_ptr<Rollback> NNIMove::applyMove(JointTree &tree) {
    auto edge = tree.getNode(nodeIndex_);
    unsigned int move = left_ ? PLL_UTREE_MOVE_NNI_LEFT : PLL_UTREE_MOVE_NNI_RIGHT;
    pll_tree_rollback_t rollback;
    assert(PLL_SUCCESS == pllmod_utree_nni(edge, move, &rollback));
    if (blo_) {
        applyBLO(tree, edge);
    }
    tree.updateBPPTree();
    return std::make_shared<Rollback>(tree, rollback);
}

ostream& NNIMove::print(ostream & os) const {
  os << "NNI(" << nodeIndex_ << ",";
  os << (left_ ? "left":"right") << ",";
  os << ")";
  return os;
}


SPRMove::SPRMove(int pruneIndex, int regraftIndex):
  pruneIndex_(pruneIndex),
  regraftIndex_(regraftIndex)
{
}

void printNode(pll_unode_s *node)
{
  cout << node;
  if (node->next) {
    cout << " " << node->next << " " << node->next->next;
  }
  cout << endl;
}

std::shared_ptr<Rollback> SPRMove::applyMove(JointTree &tree)
{
  auto prune = tree.getNode(pruneIndex_);
  auto regraft = tree.getNode(regraftIndex_);
  pll_tree_rollback_t rollback;
  assert(PLL_SUCCESS == pllmod_utree_spr(prune, regraft, &rollback));
  bool blo = true;
  if (blo) {
    applyBLO(tree, prune, regraft);
  }
  tree.updateBPPTree();
  return std::make_shared<Rollback>(tree, rollback);
}

ostream& SPRMove::print(ostream & os) const {
  os << "SPR(";
  os << pruneIndex_ << "," << regraftIndex_;
  os << ")";
  return os;
}

void SPRMove::applyBLO(JointTree &tree, 
    pll_unode_t *prune,
    pll_unode_t *regraft) 
{
    std::vector<pll_unode_t *> nodesToOptimize;
    nodesToOptimize.push_back(prune);
    if (prune->next) {
        nodesToOptimize.push_back(prune->next);
        nodesToOptimize.push_back(prune->next->next);
    }
    if (prune->back && prune->back->next) {
        nodesToOptimize.push_back(prune->back->next);
        nodesToOptimize.push_back(prune->back->next->next);
    }
    nodesToOptimize.push_back(regraft);
    if (regraft->next) {
        nodesToOptimize.push_back(regraft->next);
        nodesToOptimize.push_back(regraft->next->next);
    }
    if (regraft->back && regraft->back->next) {
        nodesToOptimize.push_back(regraft->back->next);
        nodesToOptimize.push_back(regraft->back->next->next);
    }
    optimizeBranches(tree, nodesToOptimize);
}

