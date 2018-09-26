#ifndef _JOINTTREE_H_
#define _JOINTTREE_H_

#include <Bpp/Phyl/Tree/PhyloTree.h>

#include <ale/tools/IO/IO.h>
#include <ale/tools/Utils.h>
#include <ale/tools/SpeciesGeneMapper.h>

#include <likelihoods/LibpllEvaluation.h>
#include <ale/tools/ALE/ALEevaluation.h>

#include <treeSearch/Moves.h>
#include <omp.h>

#include <sstream>
#include <stack>

using namespace std;

using BPPTree = std::shared_ptr<bpp::PhyloTree>;
using BPPNode = std::shared_ptr<bpp::PhyloNode>;
using BPPBranch = std::shared_ptr<bpp::PhyloBranch>;

int query_nni_nodes(pll_unode_t * root, vector<pll_unode_t *> &buffer);
void printLibpllNode(pll_unode_s *node, ostream &os, bool isRoot);
void printLibpllTreeRooted(pll_unode_t *root, ostream &os);
void addFromLibpll(BPPTree tree, BPPNode bppFather, pll_unode_s *libpllNode);
std::shared_ptr<bpp::PhyloTree> buildFromLibpll(std::shared_ptr<LibpllEvaluation> evaluation, pll_unode_s *libpllRoot);


class AbstractJointTree {
  public:
    virtual ~AbstractJointTree() {};
    virtual void optimizeParameters() = 0;
    virtual JointTree& getThreadInstance() = 0; 
    virtual void applyMove(shared_ptr<Move> move) = 0;
    virtual int getThreadsNumber() const { return 1; };
    virtual bool checkConsistency() {return true;}
};

class JointTree: public AbstractJointTree {
public:
    JointTree(const string &newick_file,
              const string &alignment_file,
              const string &speciestree_file,
              double dupCost,
              double lossCost):
        dupCost_(dupCost),
        lossCost_(lossCost),
        transferCost_(0.01)
    {

        info_.alignmentFilename = alignment_file;
        info_.model = "GTR";
        evaluation_ = LibpllEvaluation::buildFromFile(newick_file, info_);
        updateBPPTree();
        speciesTree_ = IO::readTreeFile(speciestree_file);
        vector<BPPTree> geneTrees(1, geneTree_);
        map_ = SpeciesGeneMapper::map(
             geneTrees.begin(), geneTrees.end(), *speciesTree_, trees)[0];
    }

    JointTree(BPPTree geneTree,
              const LibpllAlignmentInfo *alignment,
              BPPTree speciesTree,
              const SpeciesGeneMap &map,
              double dupCost,
              double lossCost):
      evaluation_(LibpllEvaluation::buildFromPhylo(geneTree, *alignment)),
      speciesTree_(speciesTree),
      map_(map),
      info_(*alignment),
      dupCost_(dupCost),
      lossCost_(lossCost),
      transferCost_(0.0)
    {
      updateBPPTree();
    }

    virtual ~JointTree() {}

    void printLibpllTree() const {
        printLibpllTreeRooted(evaluation_->getTreeInfo()->root, cout);
    }

    void printBPPTree() const {
        IO::write(*geneTree_, cout);
    }

    void printSpeciesTree() const {
        IO::write(*speciesTree_, cout);
    }

    void optimizeParameters() {
        evaluation_->optimizeAllParameters();
    }

    double computeLibpllLoglk() {
        return evaluation_->computeLikelihood();
    }

    double computeALELoglk () {
        auto genetree_copy = PhyloTreeToolBox::cloneTree(*geneTree_);
        PhyloTreeToolBox::removeArtificialGenes(*genetree_copy);
        double ale_loglk = ALEevaluation::evaluate(*genetree_copy, *speciesTree_, map_, 1, 1, dupCost_, transferCost_, lossCost_);
        return ale_loglk;
    }

    double computeJointLoglk() {
        return computeLibpllLoglk() + computeALELoglk();
    }

    void printLoglk(bool libpll = true, bool ale = true, bool joint = true, ostream &os = cout) {
        if (libpll)
            os << "libpll:\t" << computeLibpllLoglk() << " ";
        if (ale)
            os << "ale:\t\t" << computeALELoglk() << " ";
        if (joint)
            os << "joint:\t" << computeJointLoglk();
        os << endl;
    }


    // todobenoit make it faster
    pll_unode_t *getNode(int index) {
      return getTreeInfo()->subnodes[index];
    }

    void getAllNodes(vector<pll_unode_t *> &allNodes) {
        auto treeinfo = evaluation_->getTreeInfo();
        assert(query_nni_nodes(treeinfo->root, allNodes) == allNodes.size());
    }

    void getAllNodeIndices(vector<int> &allNodeIndices) {
        auto treeinfo = evaluation_->getTreeInfo();
        vector<pll_unode_t*> allNodes;
        assert(query_nni_nodes(treeinfo->root, allNodes) == allNodes.size());
        for (auto node: allNodes) 
          allNodeIndices.push_back(node->node_index);
    }

    void applyMove(shared_ptr<Move> move) {
        
        auto rollback = move->applyMove(*this);
        rollbacks_.push(rollback);
    }


    void rollbackLastMove() {
        assert(!rollbacks_.empty());
        rollbacks_.top()->applyRollback();
        rollbacks_.pop();
    }

    void save(const string &fileName) {
        ofstream os(fileName);
        IO::write(*geneTree_, os);
    }

    void updateBPPTree() {
        geneTree_ = buildFromLibpll(evaluation_, evaluation_->getTreeInfo()->root);
    }

    BPPTree getGeneTree() {
      return geneTree_;
    }

    shared_ptr<pllmod_treeinfo_t> getTreeInfo() {
        return evaluation_->getTreeInfo();
    }

    virtual JointTree& getThreadInstance() {
      return *this;
    }


private:
    std::shared_ptr<LibpllEvaluation> evaluation_;
    BPPTree geneTree_;
    BPPTree speciesTree_;
    SpeciesGeneMap map_;
    LibpllAlignmentInfo info_;
    double dupCost_;
    double lossCost_;
    double transferCost_;
    stack<shared_ptr<Rollback> > rollbacks_;
};

class ParallelJointTree: public AbstractJointTree {
public:
  ParallelJointTree(BPPTree geneTree,
    const LibpllAlignmentInfo *alignment,
    BPPTree speciesTree,
    const SpeciesGeneMap &map,
    double dupCost,
    double lossCost,
    int threads)
  {
    for (int i = 0; i < threads; ++i) {
      trees_.push_back(make_shared<JointTree>(geneTree,
            alignment, speciesTree, map, dupCost, lossCost));
    }
  }
  
  ParallelJointTree(const string &newick_file,
            const string &alignment_file,
            const string &speciestree_file,
            double dupCost,
            double lossCost,
            int threads)
  {
    for (int i = 0; i < threads; ++i) {
      trees_.push_back(make_shared<JointTree>(newick_file,
            alignment_file, speciestree_file,dupCost, lossCost));
    }
  }

  virtual ~ParallelJointTree() {}
    
  void optimizeParameters() {
    #pragma omp parallel for num_threads(getThreadsNumber())
    for (int i = 0; i < trees_.size(); ++i) {
      trees_[i]->optimizeParameters();
    }
  }

  virtual int getThreadsNumber() const {
    return trees_.size();
  }

  virtual JointTree& getThreadInstance() {
    int tid = omp_get_thread_num();
    if (tid >= trees_.size()) {
      cerr << "invalid index " << trees_.size() << " in getThreadInstance" << endl;
      exit(1);
    }
    return *trees_[tid];
  }
    
  virtual void applyMove(shared_ptr<Move> move) {
    for (auto &tree: trees_) {
      tree->applyMove(move);
    }
  }

  virtual bool checkConsistency() {
    cerr << "check consistency" << endl;
    vector<double> ll(getThreadsNumber());
    vector<vector< int>> nodesIndices(getThreadsNumber());
    #pragma omp parallel for num_threads(getThreadsNumber())
    for (int i = 0; i < trees_.size(); ++i) {
      ll[i] = trees_[i]->computeJointLoglk();
      getThreadInstance().getAllNodeIndices(nodesIndices[i]);
    }
    auto refll = ll[0];
    auto refNodesIndices = nodesIndices[0];
    for (int i = 0; i < trees_.size(); ++i) {
      if (ll[i] != refll) {
        cerr << "Error, one tree at least has a different ll" << endl;
        exit(1);
      }
      if (nodesIndices[i] != refNodesIndices) {
        cerr << "Error, one tree at least has a node indices" << endl;
        exit(1);
      }
    }
  }

  private:
    vector<shared_ptr<JointTree> > trees_;
};
#endif


