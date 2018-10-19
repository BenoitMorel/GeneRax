
#ifndef JOINTSEARCH_ALEEVALUATION_HPP
#define JOINTSEARCH_ALEEVALUATION_HPP

// Include ALE
#include <ale/tools/ALE/ALE.h>
#include <ale/tools/ALE/exODT_DL.h>
#include <ale/tools/SpeciesGeneMapper.h>
#include <ale/tools/IO/IO.h>
#include <likelihoods/ale/UndatedDLModel.hpp>

class ALEEvaluation {

public:
  ALEEvaluation(pll_rtree_t *speciesTree,
    const SpeciesGeneMap& map);

  void setRates(double dupRate, double lossRate);

  double evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo);
private:
  exODT_DL_model model;
  UndatedDLModel undatedDLModel;
};


#endif //TREERECS_ALEEVALUATION_H
