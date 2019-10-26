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
#include <IO/LibpllParsers.hpp>
#include <trees/PLLUnrootedTree.hpp>
#include <trees/PLLTreeInfo.hpp>
#include <string>
#include <memory>
#include <exception>
#include <vector>
#include <iostream>

const double TOLERANCE = 0.5;




struct LibpllAlignmentInfo {
  std::string alignmentFilename;
  std::string model;
};

/*
 * Libpll wraper to compute the phylogenetic likelihood of a tree.
 */
class LibpllEvaluation {
public:
  /**
   * Forbid copy
   */
  LibpllEvaluation(const LibpllEvaluation &) = delete;
  LibpllEvaluation & operator = (const LibpllEvaluation &) = delete;
  LibpllEvaluation(LibpllEvaluation &&) = delete;
  LibpllEvaluation & operator = (LibpllEvaluation &&) = delete;

  /*
   * Constructor 
   * @param newickStrOrFile the tree in newick format: either a string or the path to a file containing the string
   * @param isNewickAFile specifies whether newickStrOrFile is a string or a filepath
   * @param alignmentFilename path to the msa file
   * @param modelStrOrFile a std::string representing the model (GTR, DAYOFF...), or a file containing it
   */
  LibpllEvaluation(const std::string &newickStrOrFile,
      bool isNewickAFile,
      const std::string &alignmentFilename,
      const std::string &modelStrOrFile);


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
  pllmod_treeinfo_t *getTreeInfo() {return _treeInfo->getTreeInfo();}
  const pllmod_treeinfo_t *getTreeInfo() const {return _treeInfo->getTreeInfo();}

  /**
   *  Invalidate a CLV at a given node index
   *  Relevant for computeLikelihood(true) calls
   */
  void invalidateCLV(unsigned int nodeIndex);

  static void createAndSaveRandomTree(const std::string &alignmentFiilename,
    const std::string &modelStrOrFile,
    const std::string &outputTreeFile);

  std::string getModelStr();

  pll_utree_t *getGeneTree() {return _treeInfo->getTree().getRawPtr();}

private:
  /**
   * Constructors
   */
  LibpllEvaluation() {}
 

  static double optimizeAllParametersOnce(pllmod_treeinfo_t *treeinfo, double tolerance);
  
  pll_unode_t *getNode(unsigned int nodeIndex) {return _treeInfo->getTreeInfo()->subnodes[nodeIndex];}
private:
  std::unique_ptr<PLLTreeInfo> _treeInfo;
};

