// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by comten on 14/08/17.
//

#ifndef TREERECS_ALEEVALUATION_H
#define TREERECS_ALEEVALUATION_H

//#define DTL

// Include ALE
#include "ALE.h"
#ifdef DTL
#include "exODT.h"
#else
#include "exODT_DL.h"
#endif
#include "fractionMissing.h"

// Include Treerecs tools
#include <ale/tools/SpeciesGeneMapper.h>
#include <ale/tools/IO/IO.h>

/// \brief Tools used to evaluate trees with ALE.
/// \details ALE, for Amalgamated likelihood estimation by G.J. Szollosi et al. (https://github.com/ssolo/ALE),
/// is a probabilistic approach to exhaustively explore all reconciled gene trees that can be amalgamated as a
/// combination of clades observed in a sample of gene trees.
///
class ALEevaluation {
private:

protected:

public:
/// Returns ALE log likelihood of a gene tree according to a species tree and its mapping.
static double evaluate(
      const bpp::PhyloTree& genetree
      , const bpp::PhyloTree& speciestree
      , const SpeciesGeneMap& map
      , long double beta = 1
      , long double O_R = 1
      , long double delta = 0.01
      , long double tau = 0.01
      , long double lambda = 0.1
  ) {
    auto gene_tree_str = IO::PhyloTreeToNewick(genetree);
    auto species_tree_str = IO::PhyloTreeToNewick(speciestree);

    std::shared_ptr<approx_posterior> ale(new approx_posterior(gene_tree_str));

    std::vector<std::string> gene_tree_strs;
    gene_tree_strs.push_back(gene_tree_str);
    ale->observation(gene_tree_strs);

    // We initialise a coarse grained reconciliation model for calculating the sum
#ifdef DTL
    exODT_model model;
#else
    exODT_DL_model model;
#endif
    // Constructing the ALE_undated object and computing the logLk.
    model.setMap(SpeciesGeneMapper::nodeMapsToStrings(map));
    model.set_model_parameter("BOOTSTRAP_LABELS", "yes");
    model.construct_undated(species_tree_str, "");

    // Set model parameters.
    model.set_model_parameter("seq_beta", beta);
    model.set_model_parameter("O_R", O_R);
    model.set_model_parameter("delta", delta);
    model.set_model_parameter("tau", tau);
    model.set_model_parameter("lambda", lambda);

    //calculate_EGb() must always be called after changing rates to calculate E-s and G-s
    //cf. http://arxiv.org/abs/1211.4606
    model.calculate_undatedEs();
    double loglK1 = (double)std::log(model.pun(ale, false));
    double loglK = (double)std::log(model.pun(ale, false));

    return loglK;
  }

};


#endif //TREERECS_ALEEVALUATION_H
