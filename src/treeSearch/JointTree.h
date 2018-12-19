#ifndef _JOINTTREE_H_
#define _JOINTTREE_H_

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>

#include <Arguments.hpp>
#include <Logger.hpp>
#include <treeSearch/Moves.h>
#include <omp.h>

#include <parsers/GeneSpeciesMapping.hpp>

#include <sstream>
#include <stack>

using namespace std;

void printLibpllNode(pll_unode_s *node, ostream &os, bool isRoot);
void printLibpllTreeRooted(pll_unode_t *root, ostream &os);

class JointTree {
public:
    JointTree(const string &newick_string,
              const string &alignment_file,
              const string &speciestree_file,
              const string &geneSpeciesMap_file,
              Arguments::ReconciliationModel reconciliationModel,
              double dupRate,
              double lossRate);


    virtual ~JointTree() {}
    void printLibpllTree() const;
    void optimizeParameters(bool felsenstein = true, bool reconciliation = true);
    double computeLibpllLoglk(bool incremental = false);
    double computeReconciliationLoglk ();
    double computeJointLoglk();
    void printLoglk(bool libpll = true, bool rec = true, bool joint = true, Logger &os = Logger::info);
    pll_unode_t *getNode(int index);
    void applyMove(shared_ptr<Move> move);
    void optimizeMove(shared_ptr<Move> move);
  
    void invalidateCLV(pll_unode_s *node);
    void printAllNodes(ostream &os);
    void rollbackLastMove();
    void save(const string &fileName, bool append);
    shared_ptr<pllmod_treeinfo_t> getTreeInfo();
    void setRates(double dup, double loss);
    void optimizeDTRates();
    pll_rtree_t *getSpeciesTree() {return pllSpeciesTree_;}
    size_t getUnrootedTreeHash();
    shared_ptr<ReconciliationEvaluation> getReconciliationEvaluation() {return reconciliationEvaluation_;}
private:
    void findBestRates(double minDup, double maxDup,
        double minLoss, double maxLoss, int steps,
        double &bestDup,
        double &bestLoss,
        double &bestLL);
    shared_ptr<LibpllEvaluation> libpllEvaluation_;
    shared_ptr<ReconciliationEvaluation> reconciliationEvaluation_;
    pll_rtree_t *pllSpeciesTree_;
    GeneSpeciesMapping geneSpeciesMap_;
    LibpllAlignmentInfo info_;
    double dupRate_;
    double lossRate_;
    stack<shared_ptr<Rollback> > rollbacks_;
    double reconciliationLL_;
};

#endif


