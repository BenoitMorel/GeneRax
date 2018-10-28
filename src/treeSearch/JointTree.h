#ifndef _JOINTTREE_H_
#define _JOINTTREE_H_

#include <likelihoods/LibpllEvaluation.hpp>
#include <likelihoods/ALEEvaluation.hpp>

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
    JointTree(const string &newick_file,
              const string &alignment_file,
              const string &speciestree_file,
              const string &geneSpeciesMap_file,
              double dupRate,
              double lossRate);


    virtual ~JointTree() {}
    void printLibpllTree() const;
    void optimizeParameters();
    double computeLibpllLoglk();
    double computeALELoglk ();
    double computeJointLoglk();
    void printLoglk(bool libpll = true, bool ale = true, bool joint = true, Logger &os = Logger::info);
    pll_unode_t *getNode(int index);
    void applyMove(shared_ptr<Move> move);
    
    void rollbackLastMove();
    void save(const string &fileName);
    shared_ptr<pllmod_treeinfo_t> getTreeInfo();
    void setRates(double dup, double loss);
    void optimizeDTRates();
    pll_rtree_t *getSpeciesTree() {return pllSpeciesTree_;}
    int getTreeHash();
    shared_ptr<ALEEvaluation> getAleEvaluation() {return aleEvaluation_;}
private:
    shared_ptr<LibpllEvaluation> evaluation_;
    shared_ptr<ALEEvaluation> aleEvaluation_;
    pll_rtree_t *pllSpeciesTree_;
    GeneSpeciesMapping geneSpeciesMap_;
    LibpllAlignmentInfo info_;
    double dupRate_;
    double lossRate_;
    stack<shared_ptr<Rollback> > rollbacks_;
    double aleWeight_;
    double aleLL_;
};

#endif


