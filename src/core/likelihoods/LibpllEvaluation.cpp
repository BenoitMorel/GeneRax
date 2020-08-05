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


//#define RAXML_PARAM_EPSILON       0.001  //0.01
#define RAXML_BFGS_FACTOR         1e7
#define RAXML_BRLEN_SMOOTHINGS    32
#define RAXML_BRLEN_MIN           1.0e-6
#define RAXML_BRLEN_MAX           100.
#define RAXML_FREERATE_MIN        0.001
#define RAXML_FREERATE_MAX        100.
#define RAXML_BRLEN_SCALER_MIN    0.01
#define RAXML_BRLEN_SCALER_MAX    100.




std::string LibpllEvaluation::getModelStr()
{
  assign(_treeInfo->getModel(), _treeInfo->getTreeInfo()->partitions[0]);
  std::string modelStr = _treeInfo->getModel().to_string(true);
  return modelStr;
}

void LibpllEvaluation::createAndSaveRandomTree(const std::string &alignmentFilename,
    const std::string &modelStrOrFile,
    const std::string &outputTreeFile)
{
  PLLTreeInfo treeInfo("__random__", false, alignmentFilename, modelStrOrFile);
  treeInfo.getTree().save(outputTreeFile);
}
  
LibpllEvaluation::LibpllEvaluation(const std::string &newickStrOrFile,
      bool isNewickAFile,
      const std::string& alignmentFilename,
      const std::string &modelStrOrFile):
  _treeInfo(std::make_unique<PLLTreeInfo>(newickStrOrFile, isNewickAFile, alignmentFilename, modelStrOrFile))
{
}

double LibpllEvaluation::raxmlSPRRounds(unsigned int minRadius, 
    unsigned int maxRadius, 
    unsigned int thorough, 
    unsigned int toKeep, 
    double cutoff)
{
  cutoff_info_t cutoff_info;
  if (cutoff != 0.0) {
    cutoff_info.lh_dec_count = 0;
    cutoff_info.lh_dec_sum = 0.;
    cutoff_info.lh_cutoff = computeLikelihood(false) / -1000.0;
  }
  return pllmod_algo_spr_round(getTreeInfo(),
      static_cast<int>(minRadius),
      static_cast<int>(maxRadius),
      static_cast<int>(toKeep), // params.ntopol_keep
      static_cast<int>(thorough), // THOROUGH
      0, //int brlen_opt_method,
      RAXML_BRLEN_MIN,
      RAXML_BRLEN_MAX,
      RAXML_BRLEN_SMOOTHINGS,
      0.1,
      (fabs(cutoff) < std::numeric_limits<double>::epsilon()) ? 0 : &cutoff_info, //cutoff_info_t * cutoff_info,
      cutoff); //double subtree_cutoff);
}


double LibpllEvaluation::computeLikelihood(bool incremental)
{
  return pllmod_treeinfo_compute_loglh(_treeInfo->getTreeInfo(), incremental);
}

double LibpllEvaluation::optimizeBranches(double tolerance)
{
  auto toOptimize = _treeInfo->getTreeInfo()->params_to_optimize[0];
  _treeInfo->getTreeInfo()->params_to_optimize[0] = PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE;
  double res = optimizeAllParameters(tolerance);
  _treeInfo->getTreeInfo()->params_to_optimize[0] = toOptimize;
  return res;
}

double LibpllEvaluation::optimizeAllParameters(double tolerance)
{
  if (_treeInfo->getTreeInfo()->params_to_optimize[0] == 0) {
    return computeLikelihood();
  }
  double previousLogl = computeLikelihood(); 
  double newLogl = previousLogl;
  do {
    previousLogl = newLogl;
    newLogl = optimizeAllParametersOnce(_treeInfo->getTreeInfo(), tolerance);
  } while (newLogl - previousLogl > tolerance);
  return newLogl;
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
        RAXML_BRLEN_MIN,
        RAXML_BRLEN_MAX,
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
  pllmod_treeinfo_invalidate_clv(_treeInfo->getTreeInfo(), getNode(nodeIndex));
  pllmod_treeinfo_invalidate_pmatrix(_treeInfo->getTreeInfo(), getNode(nodeIndex));
}

