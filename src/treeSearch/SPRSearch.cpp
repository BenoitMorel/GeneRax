#include "SPRSearch.h"

void queryPruneIndicesRec(pll_unode_t * node,
                               vector<int> &buffer)
{
  if (node->next) {
    queryPruneIndicesRec(node->next->back, buffer);
    queryPruneIndicesRec(node->next->next->back, buffer);
    buffer.push_back(node->node_index);
  }
}

void getAllPruneIndices(JointTree &tree, vector<int> &allNodeIndices) {
  auto treeinfo = tree.getTreeInfo();
  vector<int> allNodes;
  queryPruneIndicesRec(treeinfo->root->back, allNodeIndices);
  queryPruneIndicesRec(treeinfo->root, allNodeIndices);
}


bool sprYeldsSameTree(pll_unode_t *n1, pll_unode_t *n2)
{
  return (n2 == n1) || (n2 == n1->next) || (n2 == n1->next->next)
    || (n2->back == n1) || (n2->back == n1->next) || (n2->back == n1->next->next);
}

bool isValidSPRMove(shared_ptr<pllmod_treeinfo_t> treeinfo, pll_unode_s *prune, pll_unode_s *regraft) {
  return !sprYeldsSameTree(prune, regraft);
}

bool testSPRMove(JointTree &jointTree,
    int pruneIndex,
    int regraftIndex,
    double bestLoglk,
    double &newLoglk,
    shared_ptr<Move> &newMove,
    bool check
    )
{
  if (!isValidSPRMove(jointTree.getTreeInfo(), jointTree.getNode(pruneIndex), jointTree.getNode(regraftIndex))) {
    return false;
  }
  bool res = false;
  double initialLoglk = 0.0;
  if (check) {
    initialLoglk = jointTree.computeJointLoglk();
  }
  newMove = Move::createSPRMove(pruneIndex, regraftIndex);
  jointTree.applyMove(newMove);
  newLoglk = jointTree.computeLibpllLoglk();
  if (newLoglk > bestLoglk) {
     newLoglk += jointTree.computeALELoglk();
  }
  jointTree.rollbackLastMove();
  
  if(check && fabs(initialLoglk - jointTree.computeJointLoglk()) > 0.000001) {
    cerr.precision(17);
    cerr << "rollback lead to different likelihoods: " << initialLoglk
      << " " << jointTree.computeJointLoglk() << endl;
    exit(1);
  }
  return newLoglk > bestLoglk;
}


bool SPRSearch::applySPRRound(AbstractJointTree &jointTree, int radius, double &bestLoglk) {
  cout << "SPR ROUND WITH RADIUS " << radius << endl;
  shared_ptr<Move> bestMove(0);
  vector<int> allNodes;
  getAllPruneIndices(jointTree.getThreadInstance(), allNodes);
  const size_t edgesNumber = jointTree.getThreadInstance().getTreeInfo()->tree->edge_count;
  bool foundBetterMove = false;
  
  vector<pair<int, int> > movesToExplore;
  for (int i = 0; i < allNodes.size(); ++i) {
      int pruneIndex = allNodes[i];
      vector<pll_unode_t*> regraftNodes(edgesNumber);
      unsigned int offset = 0;
      unsigned int count = 0;
      unsigned int radiusMin = 1;
      unsigned int radiusMax = radius;
      for (int r = radiusMin; r <= radiusMax; ++r) { 
          pllmod_utree_nodes_at_node_dist(jointTree.getThreadInstance().getNode(pruneIndex),
              &regraftNodes[offset],
              &count,
              r,
              r);
          offset += count;
          count = 0;
      }
      regraftNodes.resize(offset);
      vector<int> regraftIndices;
      for (auto regraftNode: regraftNodes) {
        movesToExplore.push_back(pair<int, int>(pruneIndex, regraftNode->node_index));
      }
  }


  #pragma omp parallel for num_threads(jointTree.getThreadsNumber())
  for (int i = 0; i < movesToExplore.size(); ++i) {
    if (!foundBetterMove) {
        double newLoglk;
        shared_ptr<Move> newMove;
        int pruneIndex = movesToExplore[i].first;
        int regraftIndex = movesToExplore[i].second;
        if (testSPRMove(jointTree.getThreadInstance(), pruneIndex, regraftIndex, bestLoglk, newLoglk, newMove, false)) {
        #pragma omp critical
        if (!foundBetterMove) {
          foundBetterMove = true;
          bestMove = newMove;
          bestLoglk = newLoglk;
          cout << "found a better move with loglk " << newLoglk << endl;
        }
      }
    }
  }
  if (foundBetterMove) {
    jointTree.applyMove(bestMove);
  }
  return foundBetterMove;
}


void SPRSearch::applySPRSearch(AbstractJointTree &jointTree)
{ 
  jointTree.getThreadInstance().printLoglk();
  double startingLoglk = jointTree.getThreadInstance().computeJointLoglk();
  double bestLoglk = startingLoglk;
  while (applySPRRound(jointTree, 2, bestLoglk)) {}
  while (applySPRRound(jointTree, 5, bestLoglk)) {}
  while (applySPRRound(jointTree, 10, bestLoglk)) {}
}

