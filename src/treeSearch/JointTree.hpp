#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <IO/Logger.hpp>
#include <treeSearch/Moves.hpp>
#include <IO/GeneSpeciesMapping.hpp>

#include <omp.h>
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
              const string &substitutionModel,
              const string &reconciliationModel,
              const string &reconciliationOpt,
              bool rootedGeneTree,
              bool safeMode = false,
              bool optimizeDTLRates = true,
              double dupRate = 2.0,
              double lossRate = 1.0,
              double transRate = 0.0);


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
    void printInfo();
    void rollbackLastMove();
    void save(const string &fileName, bool append);
    shared_ptr<pllmod_treeinfo_t> getTreeInfo();
    void setRates(double dup, double loss, double trans = 0);
    pll_rtree_t *getSpeciesTree() {return pllSpeciesTree_;}
    size_t getUnrootedTreeHash();
    shared_ptr<ReconciliationEvaluation> getReconciliationEvaluation() {return reconciliationEvaluation_;}
    
    pll_unode_t *getRoot() {return reconciliationEvaluation_->getRoot();}
    void setRoot(pll_unode_t * root) {reconciliationEvaluation_->setRoot(root);}
    double getDupRate() const {return dupRate_;}
    double getLossRate() const {return lossRate_;}
    double getTransferRate() const {return transRate_;}
    void inferMLScenario(Scenario &scenario) {
      reconciliationEvaluation_->inferMLScenario(libpllEvaluation_->getTreeInfo(), scenario);
    }
    bool isSafeMode() {return safeMode_;}
    void enableReconciliation(bool enable) {enableReconciliation_ = enable;}
private:
    shared_ptr<LibpllEvaluation> libpllEvaluation_;
    shared_ptr<ReconciliationEvaluation> reconciliationEvaluation_;
    pll_rtree_t *pllSpeciesTree_;
    GeneSpeciesMapping geneSpeciesMap_;
    LibpllAlignmentInfo info_;
    double dupRate_;
    double lossRate_;
    double transRate_;
    stack<shared_ptr<Rollback> > rollbacks_;
    double reconciliationLL_;
    bool optimizeDTLRates_;
    bool safeMode_;
    bool enableReconciliation_;
    string recOpt_;
};


