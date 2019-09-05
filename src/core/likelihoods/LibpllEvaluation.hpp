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
#include <iostream>

const double TOLERANCE = 0.5;


class Model;
struct pll_sequence;
using pll_sequence_ptr = std::shared_ptr<pll_sequence>;
using pll_sequences = std::vector<pll_sequence_ptr>;


struct LibpllAlignmentInfo {
  std::string alignmentFilename;
  std::string model;
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
   * @param modelStrOrFile a std::string representing the model (GTR, DAYOFF...), or a file containing it
   * @return a shared pointer wraping the LibpllEvaluation instance
   */
  static std::shared_ptr<LibpllEvaluation> buildFromString(const std::string &newickString,
      const std::string& alignmentFilename,
      const std::string &modelStrOrFile);
  static std::shared_ptr<LibpllEvaluation> buildFromFile(const std::string &newickTree,
      const LibpllAlignmentInfo &info);


  /*
   *  Compute the likelihood of the tree given the alignment
   *  @param incremental if true, only recompute invalid CLVs
   *  @return the log likelihood of the tree
   */
  double computeLikelihood(bool incremental = false);

  /**
   *  Optimize branch lengths and model parameters
   *  @return the log likeihood of the tree
   */
  double optimizeAllParameters(double tolerance = TOLERANCE);
  double optimizeBranches(double tolerance = TOLERANCE);

  double raxmlSPRRounds(int minRadius, int maxRadius, int thorough, unsigned int toKeep, double cutoff);

  /**
   *  Accessor to the wrapped treeinfo structure
   */
  std::shared_ptr<pllmod_treeinfo_t> getTreeInfo() {return _treeinfo;}

  /**
   *  Invalidate a CLV at a given node index
   *  Relevant for computeLikelihood(true) calls
   */
  void invalidateCLV(unsigned int nodeIndex);

  static void createAndSaveRandomTree(const std::string &alignmentFiilename,
    const std::string &modelStrOrFile,
    const std::string &outputTreeFile);

  std::string getModelStr();

  static bool isValidAlignment(const std::string &alignmentFilename,
      const std::string &modelStr);

  pll_utree_t *getGeneTree() {return _utree.get();}

private:
  /**
   * Constructors
   */
  LibpllEvaluation():_treeinfo(nullptr) {}
  LibpllEvaluation(const LibpllEvaluation &) = delete;
 

  /**
   * set all the null branch lenghts to length
   */
  static void setMissingBL(pll_utree_t * tree, 
    double length);

  static void parseMSA(const std::string &alignmentFilename, 
    const pll_state_t *stateMap,
    pll_sequences &sequences,
    unsigned int *&weights);

  /**
   *  parse sequences and pattern weights from fasta file
   *  @param fasta_file Input file
   *  @param stdmap state std::map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parseFasta(const char *fasta_file, 
    const pll_state_t *stateMap,
    pll_sequences &sequences,
    unsigned int *&weights);

  /**
   *  parse sequences and pattern weights from phylip file
   *  @param phylip_file Input file
   *  @param stdmap state std::map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parsePhylip(const char *phylip_file, 
    const pll_state_t *stateMap,
    pll_sequences &sequences,
    unsigned int *&weights);

  static double optimizeAllParametersOnce(pllmod_treeinfo_t *treeinfo, double tolerance);
  
  pll_unode_t *getNode(unsigned int nodeIndex) {return _treeinfo->subnodes[nodeIndex];}
private:
  std::shared_ptr<pllmod_treeinfo_t> _treeinfo;
  std::shared_ptr<pll_utree_t> _utree;
  std::shared_ptr<Model> _model;
};

