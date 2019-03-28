#pragma once

extern "C" {
#include <pll.h>
#include <pllmod_algorithm.h>
#include <pll_binary.h>
#include <pll_msa.h>
#include <pll_optimize.h>
#include <pll_tree.h>
#include <pllmod_util.h>  
}

#include <string>
#include <memory>
#include <exception>
#include <vector>

const double TOLERANCE = 0.5;

using namespace std;

struct pll_sequence;
using pll_sequence_ptr = shared_ptr<pll_sequence>;
using pll_sequences = vector<pll_sequence_ptr>;


struct LibpllAlignmentInfo {
  string alignmentFilename;
  string model;
};

/*
 * Libpll wraper to compute the phylogenetic likelihood of a tree.
 */
class LibpllEvaluation {
public:
  /*
   * Build a LibpllEvaluation instance
   * @param newickString the tree in newick format
   * @param alignmentFilename path to the msa file
   * @param modelStrOrFile a string representing the model (GTR, DAYOFF...), or a file containing it
   * @return a shared pointer wraping the LibpllEvaluation instance
   */
  static shared_ptr<LibpllEvaluation> buildFromString(const string &newickString,
      const string& alignmentFilename,
      const string &modelStrOrFile);
  static shared_ptr<LibpllEvaluation> buildFromFile(const string &newickTree,
      const LibpllAlignmentInfo &info);



  /*
   *  Compute the likelihood of the tree given the alignment
   *  @param incremental: if true, only recompute invalid CLVs
   *  @return the log likelihood of the tree
   */
  double computeLikelihood(bool incremental = false);

  /**
   *  Optimize branch lengths and model parameters
   *  @return the log likeihood of the tree
   */
  double optimizeAllParameters(double tolerance = TOLERANCE);

  double raxmlSPRRounds(int minRadius, int maxRadius, int thorough);

  /**
   *  Accessor to the wrapped treeinfo structure
   */
  shared_ptr<pllmod_treeinfo_t> getTreeInfo() {return treeinfo_;}

  /**
   *  Invalidate a CLV at a given node index
   *  Relevant for computeLikelihood(true) calls
   */
  void invalidateCLV(int nodeIndex);

  static void createAndSaveRandomTree(const string &alignmentFiilename,
    const string &modelStrOrFile,
    const string &outputTreeFile);

  string getModelStr();

private:
  /**
   * Constructors
   */
  LibpllEvaluation():treeinfo_(nullptr) {}
  LibpllEvaluation(const LibpllEvaluation &) = delete;
  
  /**
   * set all the null branch lenghts to length
   */
  static void setMissingBL(pll_utree_t * tree, 
    double length);

  /**
   *  parse sequences and pattern weights from fasta file
   *  @param fasta_file Input file
   *  @param map state map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parseFasta(const char *fasta_file, 
    const pll_state_t *map,
    pll_sequences &sequences,
    unsigned int *&weights);

  /**
   *  parse sequences and pattern weights from phylip file
   *  @param phylip_file Input file
   *  @param map state map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parsePhylip(const char *phylip_file, 
    const pll_state_t *map,
    pll_sequences &sequences,
    unsigned int *&weights);

  static double optimizeAllParametersOnce(pllmod_treeinfo_t *treeinfo, double tolerance);
  
  pll_unode_t *getNode(int nodeIndex) {return treeinfo_->subnodes[nodeIndex];}
private:
  shared_ptr<pllmod_treeinfo_t> treeinfo_;
  shared_ptr<pll_utree_t> utree_;
};

