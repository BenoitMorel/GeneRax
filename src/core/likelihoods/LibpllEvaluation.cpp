//
// Created by BenoitMorel on 23/01/18.
//

extern "C" {
  #include <pllmod_common.h>
}

#include "LibpllEvaluation.hpp"
#include <IO/LibpllParsers.hpp>
#include <map>
#include <fstream>
#include <iostream>
#include <IO/Logger.hpp>
#include <IO/Model.hpp>
#include <string>
#include <sstream>

const double DEFAULT_BL = 0.000001;
const double TOLERANCE = 0.5;


// constants taken from RAXML
#define DEF_LH_EPSILON            0.1
#define OPT_LH_EPSILON            0.1
#define RAXML_PARAM_EPSILON       0.001  //0.01
#define RAXML_BFGS_FACTOR         1e7
#define RAXML_BRLEN_SMOOTHINGS    32
#define RAXML_BRLEN_DEFAULT       0.1
#define RAXML_BRLEN_MIN           1.0e-6
#define RAXML_BRLEN_MAX           1000.
#define RAXML_BRLEN_TOLERANCE     1.0e-7
#define RAXML_FREERATE_MIN        0.001
#define RAXML_FREERATE_MAX        100.
#define RAXML_BRLEN_SCALER_MIN    0.01
#define RAXML_BRLEN_SCALER_MAX    100.

using namespace std;


struct pll_sequence {
  pll_sequence(char *label, char *seq, unsigned int len):
    label(label),
    seq(seq),
    len(len) {}
  char *label;
  char *seq;
  unsigned int len;
  ~pll_sequence() {
    free(label);
    free(seq);
  }
};

unsigned int getBestLibpllAttribute() {
  pll_hardware_probe();
  unsigned int arch = PLL_ATTRIB_ARCH_CPU;
  if (pll_hardware.avx_present) {
    arch = PLL_ATTRIB_ARCH_AVX;
  } else if (pll_hardware.sse_present) {
    arch = PLL_ATTRIB_ARCH_SSE;
  }
  arch |= PLL_ATTRIB_SITE_REPEATS;
  return  arch;
}

void utreeDestroy(pll_utree_t *utree) {
  if(!utree)
    return;
  free(utree->nodes);
  free(utree);
}

void treeinfoDestroy(pllmod_treeinfo_t *treeinfo)
{
  if (!treeinfo)
    return;
  pll_partition_destroy(treeinfo->partitions[0]);
  pll_utree_graph_destroy(treeinfo->root, 0);
  pllmod_treeinfo_destroy(treeinfo);
}

bool getNextLine(ifstream &is, string &os)
{
  while (getline(is, os)) {
    #if defined _WIN32 || defined __CYGWIN__
    os.erase(remove(os.begin(), os.end(), '\r'), os.end());
    #endif
    auto end = os.find("#");
    if (string::npos != end)
      os = os.substr(0, end);
    end = os.find(" ");
    if (string::npos != end)
      os = os.substr(0, end);
    if (os.size()) 
      return true;
  }
  return false;
}


void LibpllEvaluation::createAndSaveRandomTree(const string &alignmentFilename,
    const string &modelStrOrFile,
    const string &outputTreeFile)
{
  Logger::info << "todobenoit, fix code duplication" << endl;
  // sequences 
  pll_sequences sequences;
  unsigned int *patternWeights = nullptr;
  string modelStr = modelStrOrFile;
  ifstream f(modelStr);
  if (f.good()) {
    getline(f, modelStr);
    modelStr = modelStr.substr(0, modelStr.find(","));
  }
  Model model(modelStr);
  assert(model.num_submodels() == 1);
  pll_utree_t *utree = 0;
  if (!ifstream(alignmentFilename.c_str()).good()) {
    throw LibpllException("Alignment file " + alignmentFilename + "does not exist");
  }
  try {
    parsePhylip(alignmentFilename.c_str(),
        model.charmap(), sequences,
        patternWeights);
  } catch (...) {
    parseFasta(alignmentFilename.c_str(),
        model.charmap(), sequences, patternWeights);
  }
  vector<const char*> labels;
  for (auto seq: sequences) {
    labels.push_back(seq->label);
  }
  int seed = 0;
  utree = pllmod_utree_create_random(labels.size(), &labels[0], seed);
  LibpllParsers::saveUtree(utree->nodes[0], outputTreeFile, false);
}

shared_ptr<LibpllEvaluation> LibpllEvaluation::buildFromString(const string &newickString,
    const string& alignmentFilename,
    const string &modelStrOrFile)
{
  // sequences 
  pll_sequences sequences;
  unsigned int *patternWeights = nullptr;
  string modelStr = modelStrOrFile;
  ifstream f(modelStr);
  if (f.good()) {
    getline(f, modelStr);
    modelStr = modelStr.substr(0, modelStr.find(","));
  }
  Logger::info << "Libpll loaded model: " << modelStr << endl;
  Model model(modelStr);
  unsigned int statesNumber = model.num_states();
  assert(model.num_submodels() == 1);
  pll_utree_t *utree = 0;
  {
    if (!ifstream(alignmentFilename.c_str()).good()) {
      throw LibpllException("Alignment file " + alignmentFilename + "does not exist");
    }
    try {
      parsePhylip(alignmentFilename.c_str(),
          model.charmap(), sequences,
          patternWeights);
    } catch (...) {
      parseFasta(alignmentFilename.c_str(),
          model.charmap(), sequences, patternWeights);
    }
    // tree
    if (newickString == "__random__") {
      vector<const char*> labels;
      for (auto seq: sequences) {
        labels.push_back(seq->label);
      }
      int seed = 0;
      utree = pllmod_utree_create_random(labels.size(), &labels[0], seed);
    } else {
      utree = LibpllParsers::readNewickFromStr(newickString);
    }
  }
  // partition
  unsigned int attribute = getBestLibpllAttribute();
  unsigned int tipNumber = sequences.size();
  unsigned int innerNumber = tipNumber -1;
  unsigned int edgesNumber = 2 * tipNumber - 1;
  unsigned int sitesNumber = sequences[0]->len;
  unsigned int ratesMatrices = 1;
  pll_partition_t *partition = pll_partition_create(tipNumber,
      innerNumber,
      statesNumber,
      sitesNumber, 
      ratesMatrices, 
      edgesNumber,// prob_matrices
      model.num_ratecats(),  
      edgesNumber,// scalers
      attribute);  
  if (!partition) 
    throw LibpllException("Could not create libpll partition");
  pll_set_pattern_weights(partition, patternWeights);
  free(patternWeights);

  // fill partition
  map<string, int> tipsLabelling;
  unsigned int labelIndex = 0;
  for (auto seq: sequences) {
    tipsLabelling[seq->label] = labelIndex;
    pll_set_tip_states(partition, labelIndex, model.charmap(), seq->seq);
    labelIndex++;
  }
  sequences.clear();
  assign(partition, model);
  Logger::info << "Freq: ";
  for (int i = 0; i < 4; ++i) {
    Logger::info << partition->frequencies[0][i] << " ";
  }
  Logger::info << endl;
  pll_unode_t *root = utree->nodes[utree->tip_count + utree->inner_count - 1];
  pll_utree_reset_template_indices(root, utree->tip_count);
  setMissingBL(utree, DEFAULT_BL);
  
  // map tree to partition
  for (unsigned int i = 0; i < utree->inner_count + utree->tip_count; ++i) {
    auto node = utree->nodes[i];
    if (!node->next) { // tip!
      node->clv_index = tipsLabelling[node->label];
    }
  }
 
  // treeinfo
  int params_to_optimize = model.params_to_optimize();
  vector<unsigned int> params_indices(model.num_ratecats(), 0); 
  auto treeinfo = pllmod_treeinfo_create(root, 
      tipNumber, 1, PLLMOD_COMMON_BRLEN_SCALED);
  if (!treeinfo)
    throw LibpllException("Cannot create treeinfo");
  pllmod_treeinfo_init_partition(treeinfo, 0, partition,
      params_to_optimize,
      model.gamma_mode(),
      model.alpha(), 
      &params_indices[0],
      nullptr);
  
  shared_ptr<LibpllEvaluation> evaluation(new LibpllEvaluation());
  evaluation->treeinfo_ = shared_ptr<pllmod_treeinfo_t>(treeinfo, treeinfoDestroy); 
  evaluation->utree_ = shared_ptr<pll_utree_t>(utree, utreeDestroy); 
  return evaluation;
}
  
shared_ptr<LibpllEvaluation> LibpllEvaluation::buildFromFile(const string &newickFilename,
    const LibpllAlignmentInfo &info)
{
  ifstream t(newickFilename);
  if (!t)
    throw LibpllException("Could not load open newick file ", newickFilename);
  string str((istreambuf_iterator<char>(t)),
                       istreambuf_iterator<char>());
  return buildFromString(str,
      info.alignmentFilename,
      info.model);
}



double LibpllEvaluation::computeLikelihood(bool incremental)
{
  return pllmod_treeinfo_compute_loglh(treeinfo_.get(), incremental);
}

double LibpllEvaluation::optimizeAllParameters()
{
  Logger::timed << "Starting libpll rates optimization" << endl;
  double previousLogl = computeLikelihood(); 
  double newLogl = previousLogl;
  do {
    previousLogl = newLogl;
    newLogl = optimizeAllParametersOnce(treeinfo_.get());
  } while (newLogl - previousLogl > TOLERANCE);
  return newLogl;
}

void LibpllEvaluation::setMissingBL(pll_utree_t * tree, 
    double length)
{
  for (unsigned int i = 0; i < tree->tip_count; ++i)
    if (!tree->nodes[i]->length)
      tree->nodes[i]->length = length;
  for (unsigned int i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i) {
    if (!tree->nodes[i]->length)
      tree->nodes[i]->length = length;
    if (!tree->nodes[i]->next->length)
      tree->nodes[i]->next->length = length;
    if (!tree->nodes[i]->next->next->length)
      tree->nodes[i]->next->next->length = length;
  }  
}

void LibpllEvaluation::parseFasta(const char *fastaFile, 
    const pll_state_t *map,
    pll_sequences &sequences,
    unsigned int *&weights)
{
  auto reader = pll_fasta_open(fastaFile, pll_map_fasta);
  if (!reader) {
    throw LibpllException("Cannot parse fasta file ", fastaFile);
  }
  char * head;
  long head_len;
  char *seq;
  long seq_len;
  long seqno;
  int length;
  while (pll_fasta_getnext(reader, &head, &head_len, &seq, &seq_len, &seqno)) {
    sequences.push_back(pll_sequence_ptr(new pll_sequence(head, seq, seq_len)));
    length = seq_len;
  }
  int count = sequences.size();;
  char** buffer = (char**)malloc(count * sizeof(char *));
  for (int i = 0; i < count; ++i) {
    buffer[i] = sequences[i]->seq;
  }
  weights = pll_compress_site_patterns(buffer, map, count, &length);
  if (!weights) 
    throw LibpllException("Error while parsing fasta: cannot compress sites from ", fastaFile);
  for (int i = 0; i < count; ++i) {
    sequences[i]->len = length;
  }
  free(buffer);
  pll_fasta_close(reader);
}
  
void LibpllEvaluation::parsePhylip(const char *phylipFile, 
    const pll_state_t *map,
    pll_sequences &sequences,
    unsigned int *&weights)
{
  shared_ptr<pll_phylip_t> reader(pll_phylip_open(phylipFile, pll_map_phylip),
      pll_phylip_close);
  if (!reader) {
    throw LibpllException("Error while opening phylip file ", phylipFile);
  }
  pll_msa_t *msa = nullptr;
  // todobenoit check memory leaks when using the exception trick
  try {
    msa = pll_phylip_parse_interleaved(reader.get());
    if (!msa) {
      throw LibpllException("");
    }
  } catch (...) {
    msa = pll_phylip_parse_sequential(reader.get());
    if (!msa) {
      throw LibpllException("");
    }
  }
  weights = pll_compress_site_patterns(msa->sequence, map, msa->count, &msa->length);
  if (!weights) 
    throw LibpllException("Error while parsing fasta: cannot compress sites");
  for (auto i = 0; i < msa->count; ++i) {
    pll_sequence_ptr seq(new pll_sequence(msa->label[i], msa->sequence[i], msa->count));
    sequences.push_back(seq);
    // avoid freeing these buffers with pll_msa_destroy
    msa->label[i] = nullptr;
    msa->sequence[i] = nullptr;
  }
  pll_msa_destroy(msa);
}


double LibpllEvaluation::optimizeAllParametersOnce(pllmod_treeinfo_t *treeinfo)
{
  // This code comes from RaxML
  double new_loglh;
  unsigned int params_to_optimize = treeinfo->params_to_optimize[0]; 
  /* optimize SUBSTITUTION RATES */
  if (params_to_optimize & PLLMOD_OPT_PARAM_SUBST_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_subst_rates_treeinfo(treeinfo,
        0,
        PLLMOD_OPT_MIN_SUBST_RATE,
        PLLMOD_OPT_MAX_SUBST_RATE,
        RAXML_BFGS_FACTOR,
        RAXML_PARAM_EPSILON);
  }

  /* optimize BASE FREQS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREQUENCIES)
  {
    new_loglh = -1 * pllmod_algo_opt_frequencies_treeinfo(treeinfo,
        0,
        PLLMOD_OPT_MIN_FREQ,
        PLLMOD_OPT_MAX_FREQ,
        RAXML_BFGS_FACTOR,
        RAXML_PARAM_EPSILON);
  }

  /* optimize ALPHA */
  if (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_ALPHA,
        PLLMOD_OPT_MIN_ALPHA,
        PLLMOD_OPT_MAX_ALPHA,
        RAXML_PARAM_EPSILON);
  }

  if (params_to_optimize & PLLMOD_OPT_PARAM_PINV)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_PINV,
        PLLMOD_OPT_MIN_PINV,
        PLLMOD_OPT_MAX_PINV,
        RAXML_PARAM_EPSILON);
  }

  /* optimize FREE RATES and WEIGHTS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREE_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_rates_weights_treeinfo (treeinfo,
        RAXML_FREERATE_MIN,
        RAXML_FREERATE_MAX,
        RAXML_BFGS_FACTOR,
        RAXML_PARAM_EPSILON);

    /* normalize scalers and scale the branches accordingly */
    if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED &&
        treeinfo->partition_count > 1)
      pllmod_treeinfo_normalize_brlen_scalers(treeinfo);

  }

  if (params_to_optimize & PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE)
  {
    double brlen_smooth_factor = 0.25; // magical number from raxml
    new_loglh = -1 * pllmod_opt_optimize_branch_lengths_local_multi(treeinfo->partitions,
        treeinfo->partition_count,
        treeinfo->root,
        treeinfo->param_indices,
        treeinfo->deriv_precomp,
        treeinfo->branch_lengths,
        treeinfo->brlen_scalers,
        RAXML_BRLEN_MIN,
        RAXML_BRLEN_MAX,
        TOLERANCE,
        brlen_smooth_factor * RAXML_BRLEN_SMOOTHINGS,
        -1,  /* radius */
        1,    /* keep_update */
        PLLMOD_OPT_BLO_NEWTON_FAST,
        treeinfo->brlen_linkage,
        treeinfo->parallel_context,
        treeinfo->parallel_reduce_cb
        );
  }

  /* optimize brlen scalers, if needed */
  if (treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED &&
      treeinfo->partition_count > 1)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_BRANCH_LEN_SCALER,
        RAXML_BRLEN_SCALER_MIN,
        RAXML_BRLEN_SCALER_MAX,
        RAXML_PARAM_EPSILON);

    /* normalize scalers and scale the branches accordingly */
    pllmod_treeinfo_normalize_brlen_scalers(treeinfo);
  }
  return new_loglh;
}

void LibpllEvaluation::invalidateCLV(int nodeIndex)
{
  pllmod_treeinfo_invalidate_clv(treeinfo_.get(), getNode(nodeIndex));
  pllmod_treeinfo_invalidate_pmatrix(treeinfo_.get(), getNode(nodeIndex));
}

