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

void printLibpllNode(pll_unode_s *node, ostream &os, bool isRoot);
void printLibpllTreeRooted(pll_unode_t *root, ostream &os);
void addFromLibpll(BPPTree tree, BPPNode bppFather, pll_unode_s *libpllNode);
std::shared_ptr<bpp::PhyloTree> buildFromLibpll(std::shared_ptr<LibpllEvaluation> evaluation, pll_unode_s *libpllRoot);

class JointTree {
public:
    JointTree(const string &newick_file,
              const string &alignment_file,
              const string &speciestree_file,
              double dupRate,
              double lossRate);

    JointTree(BPPTree geneTree,
              const LibpllAlignmentInfo *alignment,
              BPPTree speciesTree,
              const SpeciesGeneMap &map,
              double dupRate,
              double lossRate);

    virtual ~JointTree() {}
    void printLibpllTree() const;
    void printBPPTree() const;
    void printSpeciesTree() const;
    void optimizeParameters();
    double computeLibpllLoglk();
    double computeALELoglk ();
    double computeJointLoglk();
    void printLoglk(bool libpll = true, bool ale = true, bool joint = true, ostream &os = cout);
    pll_unode_t *getNode(int index);
    void applyMove(shared_ptr<Move> move);
    void rollbackLastMove();
    void save(const string &fileName);
    void updateBPPTree();
    BPPTree getGeneTree();
    shared_ptr<pllmod_treeinfo_t> getTreeInfo();
    void setRates(double dup, double loss) { dupRate_ = dup; lossRate_ = loss;}
    void optimizeDTRates();
private:
    std::shared_ptr<LibpllEvaluation> evaluation_;
    BPPTree geneTree_;
    BPPTree speciesTree_;
    SpeciesGeneMap map_;
    LibpllAlignmentInfo info_;
    double dupRate_;
    double lossRate_;
    double transferRate_;
    stack<shared_ptr<Rollback> > rollbacks_;
    double aleWeight_;
};

#endif


