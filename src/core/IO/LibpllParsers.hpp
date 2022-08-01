#pragma once


#include <vector>
#include <string>
#include <IO/FamiliesFileParser.hpp>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <IO/Model.hpp>

typedef struct corax_utree_s corax_utree_t;
typedef struct corax_unode_s corax_unode_t;
typedef struct corax_rtree_s corax_rtree_t;
typedef struct corax_rnode_s corax_rnode_t;
typedef unsigned long long corax_state_t;

char * corax_rtree_export_newick(const corax_rnode_t * root,
                                   char * (*cb_serialize)(const corax_rnode_t *));

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
using SuperMatrix = std::unordered_map<std::string, std::string>;


class LibpllParsers {
public:
  LibpllParsers() = delete;

  static corax_utree_t *readNewickFromFile(const std::string &newickFile);
  static corax_utree_t *readNewickFromStr(const std::string &newickSTring);
  static corax_rtree_t *readRootedFromFile(const std::string &newickFile);
  static corax_rtree_t *readRootedFromStr(const std::string &newickFile);
  static void parseMSA(const std::string &alignmentFilename, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights);

  static unsigned int getMSALength(const std::string &alignmentFilename,
      const std::string &modelStrOrFilename);
  /**
   *  Length of MSA * ratio of non-gap characters
   */
  static double getMSAEntropy(const std::string &alignmentFilename,
      const std::string &modelStrOrFilename);
  static void fillLeavesFromUtree(corax_utree_t *utree, std::unordered_set<std::string> &leaves);
  static void fillLeavesFromRtree(corax_rtree_t *rtree, std::unordered_set<std::string> &leaves);
  
  /**
   *  return false if could not open the alignment
   */
  static bool fillLabelsFromAlignment(const std::string &alignmentFilename, 
      const std::string& modelStrOrFilename,  
      std::unordered_set<std::string> &leaves);

  static bool areLabelsValid(std::unordered_set<std::string> &leaves);
  
  static std::vector<unsigned int> parallelGetTreeSizes(const Families &families);
  static void saveUtree(const corax_unode_t *utree, 
    const std::string &fileName, 
    bool append = false);
  static void saveRtree(const corax_rnode_t *rtree, 
    const std::string &fileName);
  static void getUnodeNewickString(const corax_unode_t *rnode, std::string &newick);
  static void getRtreeHierarchicalString(const corax_rtree_t *rtree, std::string &newick);
  static std::unique_ptr<Model> getModel(const std::string &modelStrOrFilename);
  static void writeSuperMatrixFasta(const SuperMatrix &superMatrix,
      const std::string &outputFile);
private:
  /**
   *  parse sequences and pattern weights from fasta file
   *  @param fasta_file Input file
   *  @param stateMap state std::map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parseFasta(const char *fasta_file, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights);

  /**
   *  parse sequences and pattern weights from phylip file
   *  @param phylip_file Input file
   *  @param stateMap state std::map
   *  @param sequences Compressed (each site appears only once) sequences
   *  @param weights Pattern weights
   */
  static void parsePhylip(const char *phylip_file, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights);
  


};

