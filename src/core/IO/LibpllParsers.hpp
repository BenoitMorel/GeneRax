#pragma once


#include <vector>
#include <string>
#include <IO/FamiliesFileParser.hpp>
#include <unordered_set>
#include <memory>

typedef struct pll_utree_s pll_utree_t;
typedef struct pll_unode_s pll_unode_t;
typedef struct pll_rtree_s pll_rtree_t;
typedef struct pll_rnode_s pll_rnode_t;
typedef unsigned long long pll_state_t;

class Model;
struct PLLSequence {
  PLLSequence(char *label_, char *seq_, unsigned int len_):
    label(label_),
    seq(seq_),
    len(len_) {}
  char *label;
  char *seq;
  unsigned int len;
  ~PLLSequence() {
    free(label);
    free(seq);
  }
};
using PLLSequencePtr = std::unique_ptr<PLLSequence>;
using PLLSequencePtrs = std::vector<PLLSequencePtr>;

class LibpllException: public std::exception {
public:
  LibpllException(const std::string &s): msg_(s) {}
  LibpllException(const std::string &s1, 
      const std::string s2): msg_(s1 + s2) {}
  virtual const char* what() const throw() { return msg_.c_str(); }
  void append(const std::string &str) {msg_ += str;}

private:
  std::string msg_;
};

class LibpllParsers {
public:
  LibpllParsers() = delete;

  static void labelRootedTree(pll_rtree_t *tree);
  static void labelRootedTree(const std::string &unlabelledNewickFile, const std::string &labelledNewickFile);
  static pll_utree_t *readNewickFromFile(const std::string &newickFile);
  static pll_utree_t *readNewickFromStr(const std::string &newickSTring);
  static pll_rtree_t *readRootedFromFile(const std::string &newickFile);
  static pll_rtree_t *readRootedFromStr(const std::string &newickFile);
  static void parseMSA(const std::string &alignmentFilename, 
    const pll_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights);
  
  static void fillLeavesFromUtree(pll_utree_t *utree, std::unordered_set<std::string> &leaves);
  static void fillLeavesFromRtree(pll_rtree_t *rtree, std::unordered_set<std::string> &leaves);
  
  /**
   *  return false if could not open the alignment
   */
  static bool fillLabelsFromAlignment(const std::string &alignmentFilename, const std::string& modelStrOrFilename,  std::unordered_set<std::string> &leaves);

  
  static std::vector<unsigned int> parallelGetTreeSizes(const Families &families);
  static void saveUtree(const pll_unode_t *utree, 
    const std::string &fileName, 
    bool append = false);
  static void saveRtree(const pll_rnode_t *rtree, 
    const std::string &fileName);
  static void getRtreeNewickString(const pll_rtree_t *rtree, std::string &newick);
  static void getRnodeNewickString(const pll_rnode_t *rnode, std::string &newick);
  static void getRtreeHierarchicalString(const pll_rtree_t *rtree, std::string &newick);

private:
  /**
   *  parse sequences and pattern weights from fasta file
   *  @param fasta_file Input file
   *  @param stdmap state std::map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parseFasta(const char *fasta_file, 
    const pll_state_t *stateMap,
    PLLSequencePtrs &sequences,
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
    PLLSequencePtrs &sequences,
    unsigned int *&weights);
  



};

