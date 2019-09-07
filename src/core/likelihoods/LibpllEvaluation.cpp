//
// Created by BenoitMorel on 23/01/18.
//

extern "C" {
  #include <pllmod_common.h>
}

#include "LibpllEvaluation.hpp"
#include <map>
#include <fstream>
#include <iostream>
#include <IO/Logger.hpp>
#include <string>
#include <sstream>
#include <IO/Model.hpp>

const double DEFAULT_BL = 0.1;


// constants taken from RAXML
#define DEF_LH_EPSILON            0.1
#define OPT_LH_EPSILON            0.1
//#define RAXML_PARAM_EPSILON       0.001  //0.01
#define RAXML_BFGS_FACTOR         1e7
#define RAXML_BRLEN_SMOOTHINGS    32
#define RAXML_BRLEN_DEFAULT       0.1
#define RAXML_BRLEN_MIN           1.0e-6
#define RAXML_BRLEN_MAX           100.
#define RAXML_BRLEN_TOLERANCE     1.0e-7
#define RAXML_FREERATE_MIN        0.001
#define RAXML_FREERATE_MAX        100.
#define RAXML_BRLEN_SCALER_MIN    0.01
#define RAXML_BRLEN_SCALER_MAX    100.





unsigned int getBestLibpllAttribute() {
  pll_hardware_probe();
  unsigned int arch = PLL_ATTRIB_ARCH_CPU;
  if (pll_hardware.avx2_present) {
    arch = PLL_ATTRIB_ARCH_AVX2;
  } else if (pll_hardware.avx_present) {
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

bool getNextLine(std::ifstream &is, std::string &os)
{
  while (getline(is, os)) {
    #if defined _WIN32 || defined __CYGWIN__
    os.erase(remove(os.begin(), os.end(), '\r'), os.end());
    #endif
    auto end = os.find("#");
    if (std::string::npos != end)
      os = os.substr(0, end);
    end = os.find(" ");
    if (std::string::npos != end)
      os = os.substr(0, end);
    if (os.size()) 
      return true;
  }
  return false;
}

std::string LibpllEvaluation::getModelStr()
{
  assign(*_model, _treeinfo->partitions[0]);
  std::string modelStr = _model->to_string(true);
  return modelStr;
}

void LibpllEvaluation::createAndSaveRandomTree(const std::string &alignmentFilename,
    const std::string &modelStrOrFile,
    const std::string &outputTreeFile)
{
  auto evaluation = buildFromString("__random__", alignmentFilename, modelStrOrFile);
  LibpllParsers::saveUtree(evaluation->_utree->nodes[0], outputTreeFile, false);
}
  


std::shared_ptr<LibpllEvaluation> LibpllEvaluation::buildFromString(const std::string &newickString,
    const std::string& alignmentFilename,
    const std::string &modelStrOrFile)
{
  // sequences 
  PLLSequencePtrs sequences;
  unsigned int *patternWeights = nullptr;
  std::string modelStr = modelStrOrFile;
  std::ifstream f(modelStr);
  if (f.good()) {
    getline(f, modelStr);
    modelStr = modelStr.substr(0, modelStr.find(","));
  }
  Model model(modelStr);
  unsigned int statesNumber = model.num_states();
  assert(model.num_submodels() == 1);
  pll_utree_t *utree = 0;
  {
    LibpllParsers::parseMSA(alignmentFilename, model.charmap(), sequences, patternWeights);
    // tree
    if (newickString == "__random__") {
      std::vector<const char*> labels;
      for (auto seq: sequences) {
        labels.push_back(seq->label);
      }
      unsigned int seed = rand();
      utree = pllmod_utree_create_random(static_cast<unsigned int>(labels.size()), &labels[0], seed);
    } else {
      utree = LibpllParsers::readNewickFromStr(newickString);
    }
  }
  // partition
  unsigned int attribute = getBestLibpllAttribute();
  unsigned int tipNumber = static_cast<unsigned int>(sequences.size());
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
  std::map<std::string, unsigned int> tipsLabelling;
  unsigned int labelIndex = 0;
  for (auto seq: sequences) {
    tipsLabelling[seq->label] = labelIndex;
    pll_set_tip_states(partition, labelIndex, model.charmap(), seq->seq);
    labelIndex++;
  }
  sequences.clear();
  assign(partition, model);
  pll_unode_t *root = utree->nodes[utree->tip_count + utree->inner_count - 1];
  pll_utree_reset_template_indices(root, utree->tip_count);
  setMissingBL(utree, DEFAULT_BL);
  
  // std::map tree to partition
  for (unsigned int i = 0; i < utree->inner_count + utree->tip_count; ++i) {
    auto node = utree->nodes[i];
    if (!node->next) { // tip!
      node->clv_index = tipsLabelling[node->label];
    }
  }
 
  // treeinfo
  int params_to_optimize = model.params_to_optimize();
  params_to_optimize |= PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE;
  std::vector<unsigned int> params_indices(model.num_ratecats(), 0); 
  auto treeinfo = pllmod_treeinfo_create(root, 
      tipNumber, 1, PLLMOD_COMMON_BRLEN_SCALED);
  if (!treeinfo || !treeinfo->root)
    throw LibpllException("Cannot create treeinfo");
  pllmod_treeinfo_init_partition(treeinfo, 0, partition,
      params_to_optimize,
      model.gamma_mode(),
      model.alpha(), 
      &params_indices[0],
      model.submodel(0).rate_sym().data());
  
  std::shared_ptr<LibpllEvaluation> evaluation(new LibpllEvaluation());
  evaluation->_model = std::make_shared<Model>(model);
  evaluation->_treeinfo = std::shared_ptr<pllmod_treeinfo_t>(treeinfo, treeinfoDestroy); 
  evaluation->_utree = std::shared_ptr<pll_utree_t>(utree, utreeDestroy); 
  return evaluation;
}
  
std::shared_ptr<LibpllEvaluation> LibpllEvaluation::buildFromFile(const std::string &newickFilename,
    const LibpllAlignmentInfo &info)
{
  std::ifstream t(newickFilename);
  if (!t)
    throw LibpllException("Could not load open newick file ", newickFilename);
  std::string str((std::istreambuf_iterator<char>(t)),
                       std::istreambuf_iterator<char>());
  return buildFromString(str,
      info.alignmentFilename,
      info.model);
}

double LibpllEvaluation::raxmlSPRRounds(int minRadius, int maxRadius, int thorough, unsigned toKeep, double cutoff)
{
  cutoff_info_t cutoff_info;
  if (cutoff != 0.0) {
    cutoff_info.lh_dec_count = 0;
    cutoff_info.lh_dec_sum = 0.;
    cutoff_info.lh_cutoff = computeLikelihood(false) / -1000.0;
  }
  return pllmod_algo_spr_round(getTreeInfo().get(),
      minRadius,
      maxRadius,
      toKeep, // params.ntopol_keep
      thorough, // THOROUGH
      0, //int brlen_opt_method,
      RAXML_BRLEN_MIN,
      RAXML_BRLEN_MAX,
      RAXML_BRLEN_SMOOTHINGS,
      0.1,
      cutoff == 0 ? 0 : &cutoff_info, //cutoff_info_t * cutoff_info,
      cutoff); //double subtree_cutoff);
}


double LibpllEvaluation::computeLikelihood(bool incremental)
{
  return pllmod_treeinfo_compute_loglh(_treeinfo.get(), incremental);
}

double LibpllEvaluation::optimizeBranches(double tolerance)
{
  unsigned int toOptimize = _treeinfo->params_to_optimize[0];
  _treeinfo->params_to_optimize[0] = PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE;
  double res = optimizeAllParameters(tolerance);
  _treeinfo->params_to_optimize[0] = toOptimize;
  return res;
}

double LibpllEvaluation::optimizeAllParameters(double tolerance)
{
  if (_treeinfo->params_to_optimize[0] == 0) {
    return computeLikelihood();
  }
  double previousLogl = computeLikelihood(); 
  double newLogl = previousLogl;
  do {
    previousLogl = newLogl;
    newLogl = optimizeAllParametersOnce(_treeinfo.get(), tolerance);
  } while (newLogl - previousLogl > tolerance);
  return newLogl;
}

void LibpllEvaluation::setMissingBL(pll_utree_t * tree, 
    double length)
{
  for (unsigned int i = 0; i < tree->tip_count; ++i)
    if (0.0 == tree->nodes[i]->length)
      tree->nodes[i]->length = length;
  for (unsigned int i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i) {
    if (0.0 == tree->nodes[i]->length)
      tree->nodes[i]->length = length;
    if (0.0 == tree->nodes[i]->next->length)
      tree->nodes[i]->next->length = length;
    if (0.0 == tree->nodes[i]->next->next->length)
      tree->nodes[i]->next->next->length = length;
  }  
}


double LibpllEvaluation::optimizeAllParametersOnce(pllmod_treeinfo_t *treeinfo, double tolerance)
{
  // This code comes from RaxML
  double new_loglh = 0.0;
  auto params_to_optimize = treeinfo->params_to_optimize[0];
  /* optimize SUBSTITUTION RATES */
  if (params_to_optimize & PLLMOD_OPT_PARAM_SUBST_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_subst_rates_treeinfo(treeinfo,
        0,
        PLLMOD_OPT_MIN_SUBST_RATE,
        PLLMOD_OPT_MAX_SUBST_RATE,
        RAXML_BFGS_FACTOR,
        tolerance);
  }

  /* optimize BASE FREQS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREQUENCIES)
  {
    new_loglh = -1 * pllmod_algo_opt_frequencies_treeinfo(treeinfo,
        0,
        PLLMOD_OPT_MIN_FREQ,
        PLLMOD_OPT_MAX_FREQ,
        RAXML_BFGS_FACTOR,
        tolerance);
  }

  /* optimize ALPHA */
  if (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_ALPHA,
        PLLMOD_OPT_MIN_ALPHA,
        PLLMOD_OPT_MAX_ALPHA,
        tolerance);
  }

  if (params_to_optimize & PLLMOD_OPT_PARAM_PINV)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_PINV,
        PLLMOD_OPT_MIN_PINV,
        PLLMOD_OPT_MAX_PINV,
        tolerance);
  }

  /* optimize FREE RATES and WEIGHTS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREE_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_rates_weights_treeinfo (treeinfo,
        RAXML_FREERATE_MIN,
        RAXML_FREERATE_MAX,
        RAXML_BFGS_FACTOR,
        tolerance);

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
        tolerance,
        static_cast<int>(brlen_smooth_factor * static_cast<double>(RAXML_BRLEN_SMOOTHINGS)),
        -1,  /* radius */
        1,    /* keep_update */
        PLLMOD_OPT_BLO_NEWTON_SAFE,
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
        tolerance);

    /* normalize scalers and scale the branches accordingly */
    pllmod_treeinfo_normalize_brlen_scalers(treeinfo);
  }
  assert(0.0 != new_loglh);
  return new_loglh;
}

void LibpllEvaluation::invalidateCLV(unsigned int nodeIndex)
{
  pllmod_treeinfo_invalidate_clv(_treeinfo.get(), getNode(nodeIndex));
  pllmod_treeinfo_invalidate_pmatrix(_treeinfo.get(), getNode(nodeIndex));
}

